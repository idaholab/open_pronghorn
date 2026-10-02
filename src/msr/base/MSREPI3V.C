#include "MSREPI3V.h"
#include "MooseError.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace MSR
{



DenseMatrix
MSREPI3V::scaledMatrix(const DenseMatrix & A, Real scale)
{
  DenseMatrix result = A;
  for (auto & row : result)
    for (auto & value : row)
      value *= scale;
  return result;
}

std::vector<Real>
MSREPI3V::scaledVector(const std::vector<Real> & v, Real scale)
{
  std::vector<Real> result = v;
  for (auto & value : result)
    value *= scale;
  return result;
}

std::vector<Real>
MSREPI3V::add(const std::vector<Real> & a, const std::vector<Real> & b)
{
  if (a.size() != b.size())
    mooseError("MSREPI3V: vector size mismatch.");

  std::vector<Real> result(a.size(), 0.0);
  for (std::size_t i = 0; i < a.size(); ++i)
    result[i] = a[i] + b[i];
  return result;
}

std::vector<Real>
MSREPI3V::subtract(const std::vector<Real> & a, const std::vector<Real> & b)
{
  if (a.size() != b.size())
    mooseError("MSREPI3V: vector size mismatch.");

  std::vector<Real> result(a.size(), 0.0);
  for (std::size_t i = 0; i < a.size(); ++i)
    result[i] = a[i] - b[i];
  return result;
}

std::vector<Real>
MSREPI3V::matVec(const DenseMatrix & A, const std::vector<Real> & v)
{
  if (A.size() != v.size())
    mooseError("MSREPI3V: matrix/vector size mismatch.");

  std::vector<Real> result(v.size(), 0.0);
  for (std::size_t i = 0; i < A.size(); ++i)
  {
    if (A[i].size() != v.size())
      mooseError("MSREPI3V: matrix must be square.");

    for (std::size_t j = 0; j < v.size(); ++j)
      result[i] += A[i][j] * v[j];
  }
  return result;
}

bool
MSREPI3V::admissibleState(const std::vector<Real> & state,
                          bool require_nonnegative,
                          Real negativity_tolerance)
{
  if (negativity_tolerance < 0.0)
    mooseError("MSREPI3V: negativity_tolerance must be nonnegative.");

  for (const Real value : state)
    if (!std::isfinite(value) ||
        (require_nonnegative && value < -negativity_tolerance))
      return false;

  return true;
}

void
MSREPI3V::clampSmallNegatives(std::vector<Real> & state,
                              bool require_nonnegative,
                              Real negativity_tolerance)
{
  if (!require_nonnegative)
    return;

  for (auto & value : state)
    if (value < 0.0 && value >= -negativity_tolerance)
      value = 0.0;
}

Real
MSREPI3V::weightedRMSError(const std::vector<Real> & error,
                           const std::vector<Real> & y_old,
                           const std::vector<Real> & y_new,
                           Real rtol,
                           Real atol)
{
  if (error.size() != y_old.size() || error.size() != y_new.size())
    mooseError("MSREPI3V: error norm vector size mismatch.");
  if (rtol <= 0.0 || atol <= 0.0)
    mooseError("MSREPI3V: rtol and atol must be positive.");

  Real sum = 0.0;
  for (std::size_t i = 0; i < error.size(); ++i)
  {
    const Real scale = atol + rtol * std::max(std::abs(y_old[i]), std::abs(y_new[i]));
    const Real ratio = error[i] / scale;
    sum += ratio * ratio;
  }

  return std::sqrt(sum / static_cast<Real>(error.size()));
}

MSREPI3V::StepResult
MSREPI3V::attempt(const std::vector<Real> & y,
                  Real dt,
                  const RHSFunction & rhs,
                  const JacobianFunction & jacobian,
                  bool require_nonnegative,
                  Real negativity_tolerance)
{
  if (dt <= 0.0)
    mooseError("MSREPI3V: timestep must be positive.");

  const auto invalid_result = [&y]()
  {
    return StepResult{y,
                      std::vector<Real>(y.size(), std::numeric_limits<Real>::infinity()),
                      false};
  };

  if (!admissibleState(y, require_nonnegative, negativity_tolerance))
    return invalid_result();

  std::vector<Real> y_eval = y;
  clampSmallNegatives(y_eval, require_nonnegative, negativity_tolerance);

  const auto F_n = rhs(y_eval);
  const auto J_n = jacobian(y_eval);
  if (F_n.size() != y_eval.size() || J_n.size() != y_eval.size())
    mooseError("MSREPI3V: RHS/Jacobian dimension mismatch.");

  const auto hF_n = scaledVector(F_n, dt);

  const auto A_stage = scaledMatrix(J_n, 0.75 * dt);
  const auto stage_increment = MSRPhiFunctions::phi1Action(A_stage, hF_n);
  auto Y1 = add(y_eval, stage_increment);

  // Reject nonphysical stage values and clamp roundoff-scale negatives.
  if (!admissibleState(Y1, require_nonnegative, negativity_tolerance))
    return invalid_result();
  clampSmallNegatives(Y1, require_nonnegative, negativity_tolerance);

  const auto F_Y1 = rhs(Y1);
  const auto delta_Y1 = subtract(Y1, y_eval);
  const auto J_delta = matVec(J_n, delta_Y1);

  auto remainder = subtract(F_Y1, F_n);
  remainder = subtract(remainder, J_delta);

  const auto A_full = scaledMatrix(J_n, dt);
  const auto correction_input = scaledVector(remainder, 2.0 * dt);
  const auto full_actions =
      MSRPhiFunctions::phi1Phi3Actions(A_full, hF_n, correction_input);

  const auto & base_increment = full_actions.first;
  const auto & correction = full_actions.second;

  const auto y2 = add(y_eval, base_increment);
  auto y3 = add(y2, correction);

  // Reject nonphysical solutions and clamp roundoff-scale negatives.
  if (!admissibleState(y3, require_nonnegative, negativity_tolerance))
    return invalid_result();
  clampSmallNegatives(y3, require_nonnegative, negativity_tolerance);

  return {y3, correction, true};
}

std::vector<Real>
MSREPI3V::step(const std::vector<Real> & y,
               Real dt,
               const RHSFunction & rhs,
               const JacobianFunction & jacobian)
{
  return attempt(y, dt, rhs, jacobian).solution;
}

std::vector<Real>
MSREPI3V::step(const MSRChemistrySystem & chemistry,
               const std::vector<Real> & y,
               Real temperature,
               Real dt)
{
  const RHSFunction rhs = [&chemistry, temperature](const std::vector<Real> & state)
  {
    return chemistry.evaluateRHS(state, temperature);
  };

  const JacobianFunction jacobian = [&chemistry, temperature](const std::vector<Real> & state)
  {
    return chemistry.evaluateJacobian(state, temperature);
  };

  return step(y, dt, rhs, jacobian);
}

MSREPI3V::AdaptiveResult
MSREPI3V::integrateAdaptive(const std::vector<Real> & y0,
                            Real t0,
                            Real t_end,
                            Real dt_initial,
                            Real rtol,
                            Real atol,
                            const RHSFunction & rhs,
                            const JacobianFunction & jacobian,
                            bool require_nonnegative,
                            Real negativity_tolerance)
{
  if (t_end <= t0)
    mooseError("MSREPI3V: t_end must be greater than t0.");
  if (dt_initial <= 0.0)
    mooseError("MSREPI3V: initial timestep must be positive.");
  if (negativity_tolerance < 0.0)
    mooseError("MSREPI3V: negativity_tolerance must be nonnegative.");

  std::vector<Real> y = y0;
  if (!admissibleState(y, require_nonnegative, negativity_tolerance))
    mooseError("MSREPI3V: initial state is inadmissible.");
  clampSmallNegatives(y, require_nonnegative, negativity_tolerance);

  Real t = t0;
  Real dt = dt_initial;
  unsigned int accepted = 0;
  unsigned int rejected = 0;

  const unsigned int max_attempts = 100000;
  const Real safety = 0.9;
  const Real min_factor = 0.2;
  const Real max_factor = 5.0;

  for (unsigned int attempt_number = 0;
       attempt_number < max_attempts && t < t_end;
       ++attempt_number)
  {
    dt = std::min(dt, t_end - t);

    const auto trial =
        attempt(y, dt, rhs, jacobian, require_nonnegative, negativity_tolerance);

    const Real error_norm =
        trial.valid
            ? weightedRMSError(trial.error_estimate, y, trial.solution, rtol, atol)
            : std::numeric_limits<Real>::infinity();

    Real factor = max_factor;
    if (error_norm > 0.0)
      factor = safety * std::pow(error_norm, -1.0 / 3.0);
    factor = std::max(min_factor, std::min(max_factor, factor));

    if (trial.valid && error_norm <= 1.0)
    {
      y = trial.solution;
      t += dt;
      ++accepted;
      dt *= factor;
    }
    else
    {
      ++rejected;
      dt *= factor;
    }

    if (dt <= 0.0 || !std::isfinite(dt))
      mooseError("MSREPI3V: adaptive timestep became invalid.");
  }

  if (t < t_end)
    mooseError("MSREPI3V: adaptive integration exceeded the maximum number of attempts.");

  return {y, accepted, rejected, dt};
}

MSREPI3V::AdaptiveResult
MSREPI3V::integrateAdaptive(const MSRChemistrySystem & chemistry,
                            const std::vector<Real> & y0,
                            Real temperature,
                            Real t0,
                            Real t_end,
                            Real dt_initial,
                            Real rtol,
                            Real atol)
{
  const RHSFunction rhs = [&chemistry, temperature](const std::vector<Real> & state)
  {
    return chemistry.evaluateRHS(state, temperature);
  };

  const JacobianFunction jacobian = [&chemistry, temperature](const std::vector<Real> & state)
  {
    return chemistry.evaluateJacobian(state, temperature);
  };

  return integrateAdaptive(
      y0, t0, t_end, dt_initial, rtol, atol, rhs, jacobian, true, 10.0 * atol);
}

} // namespace MSR
