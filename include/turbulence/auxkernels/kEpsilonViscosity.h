//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "INSFVVelocityVariable.h"
#include "NS.h"
#include "TurbulenceMethods.h"

/**
 * Computes the turbulent viscosity for the k-Epsilon family of models.
 *
 * Supports:
 *   - Standard k-epsilon (SKE)
 *   - Standard k-epsilon Two-Layer (SKE 2L)
 *   - Standard k-epsilon Low-Re (SKE LRe)
 *   - Realizable k-epsilon (RKE)
 *   - Realizable k-epsilon Two-Layer (RKE 2L)
 *
 * Near-wall behaviour can be handled either through wall functions
 * (bulk_wall_treatment=true) or via two-layer / low-Re formulations
 * when wall distance is provided.
 */
class kEpsilonViscosity : public AuxKernel
{
public:
  static InputParameters validParams();

  void initialSetup() override;

  kEpsilonViscosity(const InputParameters & parameters);

protected:
  Real computeValue() override;

  /// The dimension of the domain
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<Real> & _u_var;
  /// y-velocity
  const Moose::Functor<Real> * _v_var;
  /// z-velocity
  const Moose::Functor<Real> * _w_var;

  /// Turbulent kinetic energy
  const Moose::Functor<Real> & _k;
  /// Turbulent kinetic energy dissipation rate
  const Moose::Functor<Real> & _epsilon;

  /// Density
  const Moose::Functor<Real> & _rho;
  /// Dynamic viscosity
  const Moose::Functor<Real> & _mu;

  /// Base C_mu closure coefficient (for standard variants and two-layer)
  const Real _C_mu;

  /// Maximum allowable mu_t_ratio : mu/mu_t
  const Real _mu_t_ratio_max;

  /// Wall boundaries
  const std::vector<BoundaryName> & _wall_boundary_names;

  /// If the user wants to enable bulk wall treatment (wall functions)
  const bool _bulk_wall_treatment;

  /// Method used for wall treatment (for bulk_wall_treatment branch)
  NS::WallTreatmentEnum _wall_treatment;

  /// Method used to limit the k-e time scale: "none" or "standard"
  const MooseEnum _scale_limiter;

  /// Whether we are using a newton solve
  const bool _newton_solve;

  // ---- Extended k-epsilon modeling options ----

  /// Selected k-epsilon model variant
  const NS::KEpsilonVariant _variant;

  /// Selected two-layer flavour (Wolfstein / Norris-Reynolds / Xu)
  const NS::TwoLayerFlavor _two_layer_flavor;

  /// Low-Re f_mu coefficients (SKE LRe)
  const Real _Cd0;
  const Real _Cd1;
  const Real _Cd2;

  /// Realizable C_mu coefficients (Ca0..Ca3)
  const Real _Ca0;
  const Real _Ca1;
  const Real _Ca2;
  const Real _Ca3;

  /// Time-scale constant C_t used in the k-epsilon time scale
  const Real _Ct;

  /// Functor giving distance to closest wall (for two-layer and low-Re)
  const Moose::Functor<Real> * _wall_distance_functor;

  /// Whether wall distance is available
  const bool _has_wall_distance;

  ///@{
  /// Maps for wall bounded elements (used only for bulk_wall_treatment)
  std::unordered_set<const Elem *> _wall_bounded;
  std::map<const Elem *, std::vector<Real>> _dist;
  std::map<const Elem *, std::vector<const FaceInfo *>> _face_infos;
  ///@}
};
