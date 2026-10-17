#include "MSRLinearImplicitMidpoint.h"

#include "libmesh/linear_implicit_system.h"

registerMooseObject("OpenPronghornApp", MSRLinearImplicitMidpoint);

InputParameters
MSRLinearImplicitMidpoint::validParams()
{
  InputParameters params = TimeIntegrator::validParams();
  params.addClassDescription(
      "Implicit midpoint written as one linear solve for linear FV systems, followed by "
      "the reconstruction y_{n+1} = 2 Y - y_n.");
  return params;
}

MSRLinearImplicitMidpoint::MSRLinearImplicitMidpoint(const InputParameters & parameters)
  : TimeIntegrator(parameters)
{
  _sys.needSolutionState(1);
}

void
MSRLinearImplicitMidpoint::computeTimeDerivatives()
{
  if (!_sys.solutionUDot())
    mooseError("MSRLinearImplicitMidpoint: Time derivative of solution (`u_dot`) is not stored.");

  NumericVector<Number> & u_dot = *_sys.solutionUDot();

  u_dot = *_solution;
  u_dot -= _solution_old;
  u_dot *= 2.0 / _dt;
  u_dot.close();

  computeDuDotDu();
}

void
MSRLinearImplicitMidpoint::computeADTimeDerivatives(ADReal & ad_u_dot,
                                                    const dof_id_type & dof,
                                                    ADReal & /*ad_u_dotdot*/) const
{
  const auto midpoint_state = ad_u_dot;
  ad_u_dot = 2.0 * (midpoint_state - _solution_old(dof)) / _dt;
}

void
MSRLinearImplicitMidpoint::postResidual(NumericVector<Number> & residual)
{
  residual += *_Re_time;
  residual += *_Re_non_time;
  residual.close();
}

Real
MSRLinearImplicitMidpoint::duDotDuCoeff() const
{
  return 2.0;
}

Real
MSRLinearImplicitMidpoint::timeDerivativeRHSContribution(
    const dof_id_type dof_id,
    const std::vector<Real> & factors) const
{
  mooseAssert(factors.size() == numStatesRequired(),
              "MSRLinearImplicitMidpoint: expected exactly one time-derivative coefficient state.");

  return 2.0 * factors[0] * _solution_old(dof_id) / _dt;
}

Real
MSRLinearImplicitMidpoint::timeDerivativeMatrixContribution(const Real factor) const
{
  return 2.0 * factor / _dt;
}

void
MSRLinearImplicitMidpoint::reconstructSolution()
{
  // [TimeIntegrator] is also instantiated on the dummy nl0 system used only
  // to satisfy MOOSE framework initialization. nl0 is not part of the
  // transport solve, so there is nothing to reconstruct there.
  if (!_linear_implicit_system)
    return;

  // The linear solve stores the midpoint stage Y in the system solution.
  // Convert it to the actual implicit-midpoint endpoint:
  //
  //   y_{n+1} = 2 Y - y_n.
  //
  // This must happen before any TIMESTEP_END transfer sends the species state
  // to the chemistry MultiApp.
  auto & solution = *_linear_implicit_system->solution;

  solution *= 2.0;
  solution.add(-1.0, _solution_old);
  solution.close();

  _linear_implicit_system->update();
}

void
MSRLinearImplicitMidpoint::postStep()
{
  // Intentionally empty.
  //
  // In MOOSE's transient execution order, nested TIMESTEP_END transfers and
  // MultiApps run inside FixedPointSolve::solveStep(), before TransientBase
  // later invokes TimeIntegrator::postStep(). Reconstructing here is therefore
  // too late for operator splitting. MSRLinearMidpointReconstruction performs
  // the reconstruction at EXEC_TIMESTEP_END before transfers instead.
}
