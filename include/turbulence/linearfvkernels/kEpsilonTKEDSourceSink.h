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
 * Kernel that adds contributions to the source and sink of the
 * turbulent dissipation rate epsilon in a linear finite-volume system.
 *
 * Responsibilities:
 *  - Bulk region:
 *      * implicit destruction term (matrix)
 *      * explicit production term (RHS), with optional model corrections
 *        and k–epsilon variant selection.
 *  - Near-wall region:
 *      * implicit Dirichlet-like treatment of epsilon using wall functions:
 *        the matrix coefficient enforces epsilon = epsilon_wall at the cell centroid.
 */
class kEpsilonTKEDSourceSink : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  kEpsilonTKEDSourceSink(const InputParameters & params);

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

  /// Turbulent kinetic energy k (TKE)
  const Moose::Functor<Real> & _k;

  /// Density
  const Moose::Functor<Real> & _rho;
  /// Dynamic viscosity
  const Moose::Functor<Real> & _mu;
  /// Turbulent viscosity mu_t
  const Moose::Functor<Real> & _mu_t;

  /// Wall boundary names
  const std::vector<BoundaryName> & _wall_boundary_names;

  /// Wall-treatment choice (equilibrium / non-equilibrium / linearized, etc.)
  const NS::WallTreatmentEnum _wall_treatment;

  /// C_mu constant used in wall functions and for length scales
  const Real _C_mu;

  // ---- epsilon model constants ----

  /// C1_epsilon (production coefficient)
  const Real _C1_eps;
  /// C2_epsilon (destruction coefficient)
  const Real _C2_eps;
  /// C3_epsilon (buoyancy coefficient)
  Real _C3_eps;
  /// C_t (used with Yap / ambient source)
  const Real _Ct;

  /// Low-Re / damping constant C in f2_SKE_LRe
  const Real _C_lowRe;
  /// Coefficient D in low-Re extra production G'
  const Real _D_lowRe;
  /// Coefficient E in low-Re extra production G'
  const Real _E_lowRe;
  /// Yap correction constant C_w
  const Real _Cw;

  // Production limiter (to match old LinearFVTKEDSourceSink behaviour)
  const Real _C_pl;

  // ---- k-epsilon variant and switches (mirrors LinearFVTKESourceSink) ----

  /// Selected k-epsilon variant
  const NS::KEpsilonVariant _variant;

  /// Switches for optional corrections
  NS::KEpsilonSwitches _switches;

  /// Turbulent Prandtl number for buoyancy (if used)
  const Real _Pr_t;

  /// Compressibility coefficient C_M (even though gamma_M mainly acts in the k-equation)
  const Real _C_M;

  /// Gravity vector (for buoyancy production)
  const libMesh::VectorValue<Real> _g;

  /// Optional functor: temperature (for buoyancy)
  const Moose::Functor<Real> * _T_functor;
  /// Optional functor: thermal expansion coefficient beta(T)
  const Moose::Functor<Real> * _beta_functor;
  /// Optional functor: speed of sound c
  const Moose::Functor<Real> * _c_functor;
  /// Optional functor: wall distance d (used for Yap and low-Re G')
  const Moose::Functor<Real> * _wall_distance_functor;

  /// Flags for optional functors
  const bool _has_T;
  const bool _has_beta;
  const bool _has_c;
  const bool _has_wall_distance;

  /// Nonlinear constitutive model
  NS::NonlinearConstitutiveRelation _nonlinear_model;

  /// Curvature/rotation correction model for Realizable variants
  NS::CurvatureCorrectionModel _curvature_model;

  ///@{
  /// Maps for wall treatment
  std::unordered_set<const Elem *> _wall_bounded;
  std::map<const Elem *, std::vector<Real>> _dist;
  std::map<const Elem *, std::vector<const FaceInfo *>> _face_infos;
  ///@}

  /**
   * Compute epsilon production "P_epsilon" in the bulk, based on the chosen
   * k–epsilon variant and switches. The returned value does *not* include
   * the C1_epsilon / T_e factor; that is applied by the kernel when assembling
   * the RHS.
   */
  Real computeBulkPe(const Moose::ElemArg & elem_arg, const Moose::StateArg & state);

  /**
   * Compute the turbulent time scale T_e used in the epsilon equation.
   * For now this is simply k / epsilon, but this function makes it easy to
   * extend to the realizable model's alternative definitions.
   */
  Real computeTimeScale(const Moose::ElemArg & elem_arg, const Moose::StateArg & state) const;

  /**
   * Compute an "epsilon at the wall" value for a wall-bounded element,
   * using simple wall function reasoning:
   *
   *   epsilon_wall ≈ C_mu^{3/4} k^{3/2} / (kappa y)
   *
   * averaged across all adjacent wall faces.
   */
  Real computeWallEpsilon(const Moose::ElemArg & elem_arg, const Moose::StateArg & state) const;
};
