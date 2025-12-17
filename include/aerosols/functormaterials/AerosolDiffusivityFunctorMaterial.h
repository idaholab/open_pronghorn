//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

/**
 * FunctorMaterial that computes Brownian and turbulent diffusion coefficients
 * for aerosols.
 *
 * Brownian diffusion (slip-corrected Stokes–Einstein):
 *
 *   D_B = k_B T C_c / (3 π μ d_p),
 *
 * where C_c is the Cunningham slip correction factor.
 *
 * Turbulent diffusion:
 *
 *   D_t = μ_t / (ρ Sc_t),
 *
 * where μ_t is the turbulent dynamic viscosity provided as a functor and Sc_t
 * is the turbulent Schmidt number.
 *
 * The material exposes three functor properties:
 *   - "D_B"      : Brownian diffusion coefficient [m^2/s]
 *   - "D_t"      : Turbulent diffusion coefficient [m^2/s]
 *   - "D_eff"    : D_B + D_t [m^2/s]
 */
class AerosolDiffusivityFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  AerosolDiffusivityFunctorMaterial(const InputParameters & parameters);

protected:
  /// Carrier fluid density [kg/m^3]
  const Moose::Functor<Real> & _rho;

  /// Carrier fluid dynamic viscosity [Pa·s]
  const Moose::Functor<Real> & _mu;

  /// Carrier fluid temperature [K]
  const Moose::Functor<Real> & _T;

  /// Turbulent dynamic viscosity [Pa·s]
  const Moose::Functor<Real> & _mu_t;

  /// Particle diameter [m]
  const Real _d_p;

  /// Gas mean free path [m]
  const Real _lambda_mfp;

  /// Turbulent Schmidt number [-]
  const Real _Sc_t;
};
