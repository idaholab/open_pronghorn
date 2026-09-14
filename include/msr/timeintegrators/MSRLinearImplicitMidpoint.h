#pragma once

#include "TimeIntegrator.h"

/**
 * Linear-system form of the implicit midpoint method.
 *
 * For a linear semi-discrete problem
 *
 *   M y' + K y = b,
 *
 * define the midpoint state Y = (y_{n+1} + y_n)/2. Then implicit midpoint is
 * equivalent to the single linear solve
 *
 *   (2 M / dt + K) Y = 2 M y_n / dt + b,
 *
 * followed by
 *
 *   y_{n+1} = 2 Y - y_n.
 *
 * LinearFVTimeDerivative obtains the 2/dt coefficients through
 * timeDerivativeMatrixContribution() and timeDerivativeRHSContribution().
 * postStep() performs the final inexpensive midpoint-state reconstruction.
 */
class MSRLinearImplicitMidpoint : public TimeIntegrator
{
public:
  static InputParameters validParams();
  MSRLinearImplicitMidpoint(const InputParameters & parameters);

  int order() override { return 2; }
  bool overridesSolve() const override { return false; }
  unsigned int numStatesRequired() const override { return 1; }

  void computeTimeDerivatives() override;
  void computeADTimeDerivatives(ADReal & ad_u_dot,
                                const dof_id_type & dof,
                                ADReal & ad_u_dotdot) const override;
  void postResidual(NumericVector<Number> & residual) override;

  /// Convert the solved midpoint stage Y into the physical end-of-step state
  /// y_{n+1} = 2 Y - y_n. This must happen before TIMESTEP_END transfers.
  void reconstructSolution();

  /// Intentionally empty: reconstruction is performed by
  /// MSRLinearMidpointReconstruction at EXEC_TIMESTEP_END, before transfers.
  void postStep() override;

  Real timeDerivativeRHSContribution(
      dof_id_type dof_id,
      const std::vector<Real> & factors = {}) const override;
  Real timeDerivativeMatrixContribution(Real factor) const override;

protected:
  Real duDotDuCoeff() const override;
};
