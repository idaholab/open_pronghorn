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
#include "NS.h"
#include "TurbulenceMethods.h"

/**
 * Functor material owning the quantities and parameters of the two-layer k-epsilon model:
 *  - the wall-proximity blending factor lambda, shared by mu_t and epsilon,
 *  - the algebraic (damped) near-wall epsilon, k^1.5 / l_eps,
 *  - the algebraic near-wall turbulent viscosity, (mu_t/mu)_2layer * mu.
 */
class kEpsilonTwoLayerFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  kEpsilonTwoLayerFunctorMaterial(const InputParameters & parameters);

protected:
  /// Turbulent kinetic energy
  const Moose::Functor<Real> & _k;
  /// Density
  const Moose::Functor<Real> & _rho;
  /// Dynamic viscosity
  const Moose::Functor<Real> & _mu;
  /// Distance to the closest wall
  const Moose::Functor<Real> & _wall_distance;

  /// Closure constant C_mu
  const Real _C_mu;
  /// Two-layer formulation
  const NS::TwoLayerFlavor _flavor;
  /// Threshold of the wall-distance Reynolds number
  const Real _Rey_star;
  /// Width of the blending region in Reynolds number
  const Real _DeltaRey;
};
