//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <vector>
#include "Moose.h"
#include "MooseUtils.h"
#include "ADReal.h"
#include "metaphysicl/raw_type.h"
#include "FEProblemBase.h"
#include "SubProblem.h"

namespace NS
{

/**
 * Enumeration of the supported k-epsilon model variants.
 *
 * These are meant to map directly from input file MooseEnum options
 * in kernels / materials.
 */
enum class KEpsilonVariant
{
  Standard,          ///< High-Reynolds-number standard k-epsilon
  StandardLowRe,     ///< Low-Reynolds-number standard k-epsilon
  StandardTwoLayer,  ///< Two-layer standard k-epsilon
  Realizable,        ///< Realizable k-epsilon
  RealizableTwoLayer ///< Two-layer realizable k-epsilon
};

/**
 * Enumeration of the two-layer model flavour.
 */
enum class TwoLayerFlavor
{
  Wolfstein,
  NorrisReynolds,
  Xu
};

/**
 * Convenience struct to hold run-time switches for optional corrections
 * that appear in a number of k-epsilon variants.
 */
struct KEpsilonSwitches
{
  bool use_buoyancy = false;             ///< Include buoyancy production G_b
  bool use_compressibility = false;      ///< Include compressibility modification gamma_M
  bool use_yap = false;                  ///< Include Yap correction gamma_y
  bool use_low_re_Gprime = false;        ///< Include low-Re extra production G'
  bool use_nonlinear = false;            ///< Include non-linear production G_nl
  bool use_curvature_correction = false; ///< Include curvature correction factor f_c
};

/**
 * Invariants of the strain-rate and rotation tensors together with the velocity divergence.
 *
 * Definitions:
 *   S_ij = 0.5 (du_i/dx_j + du_j/dx_i),
 *   W_ij = 0.5 (du_i/dx_j - du_j/dx_i),
 *   S2   = 2 S_ij S_ij,
 *   W2   = 2 W_ij W_ij.
 *
 * The factor of 2 in S2 and W2 is convenient because the magnitude often
 * appears as |S| = sqrt(2 S_ij S_ij).
 */
struct StrainRotationInvariants
{
  Real S2;    ///< 2 S_ij S_ij
  Real W2;    ///< 2 W_ij W_ij
  Real div_u; ///< divergence of velocity
};

/**
 * Compute the strain / rotation invariants and divergence from the velocity functors.
 *
 * This function is dimension-agnostic:
 *  - v == nullptr implies 1D,
 *  - w == nullptr and v != nullptr implies 2D,
 *  - v != nullptr and w != nullptr implies 3D.
 */
StrainRotationInvariants
computeStrainRotationInvariants(const Moose::Functor<Real> & u,
                                const Moose::Functor<Real> * v,
                                const Moose::Functor<Real> * w,
                                const Moose::ElemArg & elem_arg,
                                const Moose::StateArg & state);

/**
 * Small helper to compute the two-layer constant c_l from C_mu.
 * c_l = 0.42 * C_mu^{-3/4}
 */
Real cl_from_Cmu(const Real Cmu);

/**
 * Two-layer length scales (epsilon length and viscosity ratio).
 */
struct TwoLayerLengths
{
  Real l_eps;    ///< length scale used for epsilon, l_epsilon
  Real mu_ratio; ///< (mu_t / mu)_2layer
};

/**
 * Wolfstein two-layer model for l_eps and (mu_t/mu)_2layer.
 */
TwoLayerLengths
twoLayerWolfstein(const Real Cmu, const Real d, const Real Re_d);

/**
 * Norris–Reynolds two-layer model for l_eps and (mu_t/mu)_2layer.
 */
TwoLayerLengths
twoLayerNorrisReynolds(const Real Cmu, const Real d, const Real Re_d);

/**
 * Xu natural-convection two-layer model variant.
 *
 * yv_star is the dimensionless wall distance used in the Xu model.
 * The expressions follow the STAR-CCM+ documentation closely.
 */
TwoLayerLengths
twoLayerXu(const Real Cmu, const Real d, const Real Re_d, const Real yv_star);

/**
 * Low-Re standard k-epsilon damping function f2.
 *
 * f2 = 1 - C exp(-Re_t^2)
 */
Real
f2_SKE_LRe(const Real C, const Real Re_t);

/**
 * Low-Re standard k-epsilon damping function for the turbulent viscosity.
 *
 * fmu = 1 - exp( - (Cd0 sqrt(Re_d) + Cd1 Re_d + Cd2 Re_d^2) )
 */
Real
fmu_SKE_LRe(const Real Cd0, const Real Cd1, const Real Cd2, const Real Re_d);

/**
 * Realizable k-epsilon f2 function.
 *
 * f2 = k / (k + sqrt(nu epsilon))
 */
Real
f2_RKE(const Real k, const Real nu, const Real eps);

/**
 * Realizable k-epsilon C_mu expression based on invariants of S and W.
 *
 * Cmu = Ca0 / (Ca1 + Ca2 S_bar + Ca3 W_bar)
 * where S_bar = (k/eps) sqrt(S2) and W_bar = (k/eps) sqrt(W2)
 * with S2 = 2 S_ij S_ij and W2 = 2 W_ij W_ij.
 */
Real
Cmu_realizable(const Real Ca0,
               const Real Ca1,
               const Real Ca2,
               const Real Ca3,
               const Real S2,
               const Real W2,
               const Real k,
               const Real eps);

/**
 * Turbulent production G_k with optional compressibility terms.
 *
 * Gk = mu_t S^2
 *    - 2/3 rho k div(u) (if include_compressibility_terms)
 *    - 2/3 mu_t (div(u))^2 (if include_compressibility_terms)
 */
Real
computeGk(const Real mu_t,
          const Real S2,
          const Real rho,
          const Real k,
          const Real div_u,
          const bool include_compressibility_terms);

/**
 * Buoyancy production G_b.
 *
 * Gb = beta mu_t / Pr_t * (gradT · g)
 */
Real
computeGb(const Real beta,
          const Real mu_t,
          const Real Pr_t,
          const libMesh::VectorValue<Real> & grad_T,
          const libMesh::VectorValue<Real> & g);

/**
 * Compressibility modification gamma_M.
 *
 * gamma_M = rho C_M k epsilon / c^2
 */
Real
computeGammaM(const Real rho,
              const Real C_M,
              const Real k,
              const Real eps,
              const Real c);

/**
 * Yap correction gamma_y.
 *
 * gamma_y = C_w epsilon^2 / k * max[(l / l_eps - 1) (l / l_eps)^2, 0]
 */
Real
computeGammaY(const Real C_w,
              const Real eps,
              const Real k,
              const Real l,
              const Real l_eps);

/**
 * Low-Re additional production G' for the standard k-epsilon model.
 *
 * G' = D f2 (Gk + 2 mu_t k / d^2) exp(-E Re_d^2)
 */
Real
computeGprime(const Real D,
              const Real E,
              const Real f2,
              const Real Gk,
              const Real mu_t,
              const Real k,
              const Real d,
              const Real Re_d);

} // namespace NS