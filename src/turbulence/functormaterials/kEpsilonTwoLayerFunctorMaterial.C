//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "kEpsilonTwoLayerFunctorMaterial.h"

registerMooseObject("OpenPronghornApp", kEpsilonTwoLayerFunctorMaterial);

InputParameters
kEpsilonTwoLayerFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Computes the quantities of the two-layer k-epsilon model: the wall-proximity blending "
      "factor lambda, the algebraic near-wall epsilon and the algebraic near-wall turbulent "
      "viscosity.");

  params.addRequiredParam<MooseFunctorName>(NS::TKE, "Coupled turbulent kinetic energy.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Fluid density.");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("wall_distance", "Distance to the closest wall.");

  params.addParam<Real>("C_mu", 0.09, "Turbulent kinetic energy closure constant C_mu.");

  MooseEnum two_layer_flavor("Wolfstein NorrisReynolds Xu", "Wolfstein");
  params.addParam<MooseEnum>("two_layer_flavor",
                             two_layer_flavor,
                             "Two-layer formulation to use (Wolfstein, NorrisReynolds, Xu).");

  params.addParam<Real>("Rey_star",
                        60.0,
                        "Wall-distance Reynolds number below which the flow is considered "
                        "viscosity affected.");
  params.addParam<Real>("DeltaRey", 10.0, "Width of the blending region in Reynolds number.");

  params.addParam<MooseFunctorName>(
      "two_layer_lambda", "two_layer_lambda", "Name of the blending factor functor (0..1).");
  params.addParam<MooseFunctorName>(
      "two_layer_epsilon", "two_layer_epsilon", "Name of the algebraic near-wall epsilon functor.");
  params.addParam<MooseFunctorName>("two_layer_mu_t",
                                    "two_layer_mu_t",
                                    "Name of the algebraic near-wall turbulent viscosity functor.");
  return params;
}

kEpsilonTwoLayerFunctorMaterial::kEpsilonTwoLayerFunctorMaterial(const InputParameters & params)
  : FunctorMaterial(params),
    _k(getFunctor<Real>(NS::TKE)),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _wall_distance(getFunctor<Real>("wall_distance")),
    _C_mu(getParam<Real>("C_mu")),
    _flavor(getParam<MooseEnum>("two_layer_flavor").getEnum<NS::TwoLayerFlavor>()),
    _Rey_star(getParam<Real>("Rey_star")),
    _DeltaRey(getParam<Real>("DeltaRey"))
{
  // Wall-distance Reynolds number and two-layer lengths, shared by the three functors
  auto lengths = [this](const auto & r, const auto & t)
  {
    const Real k = std::max(_k(r, t), 0.0);
    const Real nu = _mu(r, t) / _rho(r, t);
    const Real d = _wall_distance(r, t);
    const Real Red = std::sqrt(k) * d / std::max(nu, 1e-12);
    return std::make_pair(Red, NS::twoLayerLengths(_flavor, _C_mu, d, Red, nu, k));
  };

  addFunctorProperty<Real>(
      "two_layer_lambda",
      [this, lengths](const auto & r, const auto & t)
      { return NS::twoLayerBlendingFactor(lengths(r, t).first, _Rey_star, _DeltaRey); });

  addFunctorProperty<Real>("two_layer_epsilon",
                           [this, lengths](const auto & r, const auto & t)
                           {
                             const Real k = std::max(_k(r, t), 0.0);
                             const auto tl = lengths(r, t).second;
                             return std::pow(k, 1.5) / std::max(tl.l_eps, 1e-20);
                           });

  addFunctorProperty<Real>("two_layer_mu_t",
                           [this, lengths](const auto & r, const auto & t)
                           { return lengths(r, t).second.mu_ratio * _mu(r, t); });
}
