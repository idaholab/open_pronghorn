#pragma once

#include "MSRChemistrySystem.h"
#include "MSRPhiFunctions.h"

#include <functional>
#include <vector>

namespace MSR
{

class MSREPI3V
{
public:
  using RHSFunction = std::function<std::vector<Real>(const std::vector<Real> &)>;
  using JacobianFunction = std::function<DenseMatrix(const std::vector<Real> &)>;

  struct StepResult
  {
    std::vector<Real> solution;
    std::vector<Real> error_estimate;
    bool valid = true;
  };

  struct AdaptiveResult
  {
    std::vector<Real> solution;
    unsigned int accepted_steps = 0;
    unsigned int rejected_steps = 0;
    Real final_dt = 0.0;
  };

  static std::vector<Real> step(const MSRChemistrySystem & chemistry,
                                const std::vector<Real> & y,
                                Real temperature,
                                Real dt);

  static std::vector<Real> step(const std::vector<Real> & y,
                                Real dt,
                                const RHSFunction & rhs,
                                const JacobianFunction & jacobian);

  /// Perform one EPI3V step. Small negative values within negativity_tolerance
  /// are clamped to zero when nonnegative states are required.
  static StepResult attempt(const std::vector<Real> & y,
                            Real dt,
                            const RHSFunction & rhs,
                            const JacobianFunction & jacobian,
                            bool require_nonnegative = false,
                            Real negativity_tolerance = 0.0);

  static AdaptiveResult integrateAdaptive(const std::vector<Real> & y0,
                                          Real t0,
                                          Real t_end,
                                          Real dt_initial,
                                          Real rtol,
                                          Real atol,
                                          const RHSFunction & rhs,
                                          const JacobianFunction & jacobian,
                                          bool require_nonnegative = false,
                                          Real negativity_tolerance = 0.0);

  static AdaptiveResult integrateAdaptive(const MSRChemistrySystem & chemistry,
                                          const std::vector<Real> & y0,
                                          Real temperature,
                                          Real t0,
                                          Real t_end,
                                          Real dt_initial,
                                          Real rtol,
                                          Real atol);

private:
  static DenseMatrix scaledMatrix(const DenseMatrix & A, Real scale);
  static std::vector<Real> scaledVector(const std::vector<Real> & v, Real scale);
  static std::vector<Real> add(const std::vector<Real> & a,
                               const std::vector<Real> & b);
  static std::vector<Real> subtract(const std::vector<Real> & a,
                                    const std::vector<Real> & b);
  static std::vector<Real> matVec(const DenseMatrix & A,
                                  const std::vector<Real> & v);
  static Real weightedRMSError(const std::vector<Real> & error,
                               const std::vector<Real> & y_old,
                               const std::vector<Real> & y_new,
                               Real rtol,
                               Real atol);
  static bool admissibleState(const std::vector<Real> & state,
                              bool require_nonnegative,
                              Real negativity_tolerance);
  static void clampSmallNegatives(std::vector<Real> & state,
                                  bool require_nonnegative,
                                  Real negativity_tolerance);
};

} // namespace MSR
