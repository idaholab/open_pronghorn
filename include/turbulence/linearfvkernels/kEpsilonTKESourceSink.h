//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"
#include "NS.h"
#include "TurbulenceMethods.h"

#include <unordered_set>
#include <map>
#include <vector>

/**
 * Kernel that adds contributions to the source and sink of turbulent kinetic
 * energy (TKE) in a linear finite-volume system.
 *
 * Responsibilities:
 *  - Bulk region:
 *      * implicit destruction term (matrix)
 *      * explicit production term (RHS), with optional model corrections
 *  - Near-wall region:
 *      * implicit production and destruction using wall functions (matrix)
 *      * no explicit RHS contribution
 *
 * This class extends the standard MOOSE Navier-Stokes implementation by
 * allowing k-epsilon model variants and optional corrections to be configured
 * through input parameters, but the near-wall formulation is shared across
 * variants for now.
 */
class kEpsilonTKESourceSink : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  /// Constructor
  kEpsilonTKESourceSink(const InputParameters & params);

  /// Setup tasks that require access to the mesh/problem
  virtual void initialSetup() override;

  /// Volumetric contribution to the system matrix
  virtual Real computeMatrixContribution() override;

  /// Volumetric contribution to the right-hand side
  virtual Real computeRightHandSideContribution() override;

protected:
  /// Convenience: domain dimension
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<Real> & _u_var;
  /// y-velocity (may be null in 1D)
  const Moose::Functor<Real> * _v_var;
  /// z-velocity (may be null in 1D/2D)
  const Moose::Functor<Real> * _w_var;

  /// Dissipation rate epsilon (TKED)
  const Moose::Functor<Real> & _epsilon;

  /// Density
  const Moose::Functor<Real> & _rho;
  /// Dynamic viscosity
  const Moose::Functor<Real> & _mu;
  /// Turbulent viscosity μ_t
  const Moose::Functor<Real> & _mu_t;

  /// Wall boundary names
  const std::vector<BoundaryName> & _wall_boundary_names;

  /// Wall-treatment choice (equilibrium / non-equilibrium / linearized, etc.)
  const NS::WallTreatmentEnum _wall_treatment;

  /// C_μ constant used in the wall functions
  const Real _C_mu;

  /// Production limiter constant multiplier
  const Real _C_pl;

  // ---- Extension points for turbulence modeling (bulk production) ----

  /// Selected k-epsilon variant (standard, realizable, two-layer, etc.)
  const NS::KEpsilonVariant _variant;

  /// Switches controlling optional corrections in the production term
  NS::KEpsilonSwitches _switches;

  /// Turbulent Prandtl number for buoyancy (if used)
  const Real _Pr_t;

  /// Constant used in compressibility correction γ_M (if used)
  const Real _C_M;

  /// Gravity vector (for buoyancy production), if used
  const libMesh::VectorValue<Real> _g;

  /// Optional functor: temperature (for buoyancy)
  const Moose::Functor<Real> * _T_functor;
  /// Optional functor: thermal expansion coefficient beta(T)
  const Moose::Functor<Real> * _beta_functor;
  /// Optional functor: speed of sound c
  const Moose::Functor<Real> * _c_functor;

  /// Flags indicating whether the optional fields above are actually present
  const bool _has_T;
  const bool _has_beta;
  const bool _has_c;

  /// Nonlinear constitutive model
  NS::NonlinearConstitutiveRelation _nonlinear_model;

  /// Curvature/rotation correction model for Realizable variants
  NS::CurvatureCorrectionModel _curvature_model;

  /// Per-element flag: is this element adjacent to at least one wall?
  std::unordered_set<const Elem *> _wall_bounded;

  /// For each wall-bounded element: distances from centroid to the wall faces
  std::map<const Elem *, std::vector<Real>> _dist;

  /// For each wall-bounded element: pointers to the wall faces
  std::map<const Elem *, std::vector<const FaceInfo *>> _face_infos;

  /**
   * Compute the bulk production term G_k for the current element,
   * including selected turbulence-model corrections.
   *
   * This function is only used for elements that are *not* in the near-wall
   * set; near-wall production/destruction is handled in computeWallTerms.
   */
  Real computeBulkProduction(const Moose::ElemArg & elem_arg, const Moose::StateArg & state) const;

  /**
   * Compute near-wall production and destruction coefficients for the current
   * wall-bounded element, returned as {production, destruction}.
   *
   * Keeping them separate guarantees a positive matrix diagonal (destruction only
   * in the matrix) while moving the bounded production to the explicit RHS.
   * Previously mixing them as (destruction − production) in the matrix caused a
   * negative diagonal whenever wall production exceeded destruction (large velocity
   * gradients / poor mesh quality), which was the primary driver of TKE → ∞.
   */
  std::pair<Real, Real> computeWallTerms(const Moose::ElemArg & elem_arg,
                                         const Moose::StateArg & state) const;

  // ---- Robustness controls ----

  /// Minimum k for destruction coefficient guard (ρε/max(k,k_min)).
  const Real _k_min;

  /// Maximum ε read from the TKED functor when computing the bulk destruction ρε/k.
  /// Prevents large algebraic near-wall ε values from destroying k in adjacent bulk cells.
  const Real _eps_functor_max;

  /// Minimum ε used inside the C_pl production limiter.
  /// Prevents the limiter from being trivially zero when ε ≈ 0 at initialisation.
  const Real _eps_min;

  /// Maximum μ_t/μ ratio applied inside the production and wall-shear computations.
  /// Guards against stale over-large μ_t values when k is large but ε has not yet
  /// caught up (common in the first few SIMPLE outer iterations).
  const Real _mu_t_prod_max;

  /// Durbin realizability coefficient C_pk (default 0 = disabled).
  /// When > 0, production is additionally bounded by  C_pk · ρ · k · |S|,
  /// a k-based limit that is effective even when ε ≈ 0 at start-up.
  /// Recommended: 0.667 (= 2/3).
  const Real _C_pk;

  /// Physical lower bound on TKE. If the cell k drops below this value a strong
  /// penalty source overrides the normal physics to drive k back above the floor.
  const Real _tke_min_phys;

  /// Physical upper bound on TKE. Analogous ceiling enforcement.
  const Real _tke_max_phys;

  /// If true, use the Kato–Launder (1993) production form G_k = μ_t |S| |Ω|.
  const bool _use_kato_launder;

  /// Velocity gradient method for the turbulence production term.
  NS::TurbVelocityGradientMethod _grad_method;

  /// Penalty coefficient used for physical bounds enforcement.
  /// Must be large enough to dominate the normal matrix diagonal.
  static constexpr Real _bounds_penalty = 1e8;
};
