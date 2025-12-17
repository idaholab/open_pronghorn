//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVAdvectionDiffusionFunctorRobinBC.h"

/**
 * Linear FV boundary condition that applies aerosol deposition at walls using
 * the Lai & Nazaroff (2000) drift-flux-based deposition velocity model.
 *
 * The BC is expressed in Robin form
 *
 *   alpha * grad(C)·n + beta * C = gamma
 *
 * with
 *
 *   alpha = 1
 *   beta  = V_d / D
 *   gamma = 0
 *
 * so that the diffusive flux at the wall is
 *
 *   J_dep = -D * grad(C)·n = V_d * C
 *
 * where C is the transported scalar (e.g. aerosol mass or volume fraction)
 * and V_d is the deposition velocity, which depends on:
 *   - friction velocity u_*,
 *   - Brownian diffusivity (through the Schmidt number),
 *   - gravitational settling velocity normal to the wall,
 *   - wall orientation (floor, ceiling, vertical).
 */
class LinearFVAerosolDepositionBC : public LinearFVAdvectionDiffusionFunctorRobinBC
{
public:
  static InputParameters validParams();

  LinearFVAerosolDepositionBC(const InputParameters & parameters);

protected:
  /// Robin coefficients alpha, beta, gamma
  virtual Real getAlpha(Moose::FaceArg face, Moose::StateArg state) const override;
  virtual Real getBeta(Moose::FaceArg face, Moose::StateArg state) const override;
  virtual Real getGamma(Moose::FaceArg face, Moose::StateArg state) const override;

  /// Compute the deposition velocity for the current face [m/s]
  Real computeDepositionVelocity(const Moose::FaceArg & face) const;

  /// Ceiling (upward facing) deposition velocity (Lai & Nazaroff)
  Real computeCeilingDepositionVelocity(Real Vgs, Real u_star) const;

  /// Floor (downward facing) deposition velocity (Lai & Nazaroff)
  Real computeFloorDepositionVelocity(Real Vgs, Real u_star) const;

  /// Vertical-wall deposition velocity (Lai & Nazaroff)
  Real computeVerticalDepositionVelocity(Real u_star, Real mu, Real rho, Real T) const;

  /// Integral I in Lai & Nazaroff model (depends on Sc and r+)
  Real computeIntegralI(Real u_star, Real mu, Real rho, Real T) const;

  /// Cunningham slip correction factor (same correlation as in LinearFVAerosolDriftFlux)
  Real computeCunninghamSlipCorrection() const;

  /// Brownian diffusivity [m^2/s] using the slip-corrected Stokes–Einstein relation
  Real computeBrownianDiffusivity(Real mu, Real T) const;

  /// Gravitational settling velocity vector at the face [m/s]
  RealVectorValue computeGravitationalSettlingVelocity(Real mu, Real rho) const;

  // ---- Functors / material properties ----

  /// Carrier fluid density [kg/m^3]
  const Moose::Functor<Real> & _rho;

  /// Carrier fluid dynamic viscosity [Pa·s]
  const Moose::Functor<Real> & _mu;

  /// Carrier fluid temperature [K]
  const Moose::Functor<Real> & _T;

  /// Friction velocity u_* at the wall [m/s]
  const Moose::Functor<Real> & _u_star;

  // ---- Particle properties ----

  /// Particle diameter [m]
  const Real _d_p;

  /// Particle density [kg/m^3]
  const Real _rho_p;

  /// Gas mean free path [m]
  const Real _lambda_mfp;

  /// Gravitational acceleration vector [m/s^2]
  const RealVectorValue _gravity;

  /// The functor for the diffusion coefficient (must match the diffusion used in the PDE)
  const Moose::Functor<Real> & _diffusion_coeff;
};
