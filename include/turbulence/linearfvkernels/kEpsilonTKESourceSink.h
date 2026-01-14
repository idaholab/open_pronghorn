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
   * set; near-wall production/destruction is handled in computeMatrixContribution.
   */
  Real computeBulkProduction(const Moose::ElemArg & elem_arg, const Moose::StateArg & state) const;
};
