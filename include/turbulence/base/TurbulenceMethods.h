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
#include "RankTwoTensor.h"
#include "libmesh/elem.h"

namespace NS
{

/**
 * NOTE ON GRADIENT QUALITY IN TURBULENCE MODELS
 * -----------------------------------------------
 * All production terms (G_k, S2, W2, ...) are computed from the velocity gradient
 * tensor obtained via MOOSE's functor `.gradient()` interface. The underlying
 * gradient reconstruction method (Green-Gauss or Least-Squares) is configured at
 * the variable level via the mesh/system settings and is shared with the
 * Navier-Stokes momentum equations.
 *
 * For improved robustness in cells with ill-conditioned geometry or near
 * stagnation points, it is recommended to select the Least-Squares gradient
 * method (or a Mixed/Gradient-Limited scheme) for the velocity variables in
 * the input file.  Locally overriding the gradient method inside these kernels
 * is not practical without bypassing the MOOSE functor framework entirely, which
 * would require custom neighbor-cell iteration and is intentionally avoided here.
 *
 * The S2 and W2 invariants are clamped to [0, ∞) at every evaluation point to
 * prevent NaN propagation from floating-point noise in ill-conditioned cells.
 * Production is further bounded by the C_pl limiter (C_pl * rho * eps) to prevent
 * blow-up during early iterations when the background solution is far from converged.
 */

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
 * Enumeration of the nonlinear closure options.
 */
enum class NonlinearConstitutiveRelation
{
  None,
  Quadratic,
  Cubic
};

enum class CurvatureCorrectionModel
{
  None,
  Standard // Spalart–Shur-style rotation/curvature correction
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
StrainRotationInvariants computeStrainRotationInvariants(const Moose::Functor<Real> & u,
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
 *
 * The viscosity ratio is defined as:
 *   mu_ratio = mu_t^{2L} / mu  = C_mu^{1/4} * c_l * Re_d * f_mu(Re_d)
 *
 * where c_l = 0.42 * C_mu^{-3/4} and f_mu is a variant-specific damping function.
 * This matches STAR-CCM+ (Eq. 1074 Wolfstein / 1076 Norris-Reynolds).
 *
 * Note: the (mu_t / mu) ratio scales linearly with Re_d (first power), not Re_d^{1/4}.
 */
struct TwoLayerLengths
{
  Real l_eps;    ///< length scale for epsilon: l_eps = c_l * d * f_eps(Re_d)
  Real mu_ratio; ///< (mu_t / mu)_2layer = C_mu^{1/4} * c_l * Re_d * f_mu(Re_d)
};

/**
 * Wolfstein two-layer model for l_eps and (mu_t/mu)_2layer.
 */
TwoLayerLengths twoLayerWolfstein(const Real Cmu, const Real d, const Real Re_d);

/**
 * Norris–Reynolds two-layer model for l_eps and (mu_t/mu)_2layer.
 */
TwoLayerLengths twoLayerNorrisReynolds(const Real Cmu, const Real d, const Real Re_d);

/**
 * Xu natural-convection two-layer model variant.
 *
 * yv_star is the dimensionless wall distance used in the Xu model.
 * The expressions follow the STAR-CCM+ documentation closely.
 */
TwoLayerLengths twoLayerXu(const Real Cmu, const Real d, const Real Re_d, const Real yv_star);

/**
 * Low-Re standard k-epsilon damping function f2.
 *
 * f2 = 1 - C exp(-Re_t^2)
 */
Real f2_SKE_LRe(const Real C, const Real Re_t);

/**
 * Low-Re standard k-epsilon damping function for the turbulent viscosity.
 *
 * fmu = 1 - exp( - (Cd0 sqrt(Re_d) + Cd1 Re_d + Cd2 Re_d^2) )
 */
Real fmu_SKE_LRe(const Real Cd0, const Real Cd1, const Real Cd2, const Real Re_d);

/**
 * Realizable k-epsilon f2 function.
 *
 * f2 = k / (k + sqrt(nu epsilon))
 */
Real f2_RKE(const Real k, const Real nu, const Real eps);

/**
 * Realizable k-epsilon C_mu expression based on invariants of S and W.
 *
 * Cmu = Ca0 / (Ca1 + Ca2 S_bar + Ca3 W_bar)
 * where S_bar = (k/eps) sqrt(S2) and W_bar = (k/eps) sqrt(W2)
 * with S2 = 2 S_ij S_ij and W2 = 2 W_ij W_ij.
 */
Real Cmu_realizable(const Real Ca0,
                    const Real Ca1,
                    const Real Ca2,
                    const Real Ca3,
                    const Real S2,
                    const Real W2,
                    const Real k,
                    const Real eps);

/**
 * Compute the variable C_mu used with non-linear constitutive relations.
 *
 *   C_mu = Ca0 / (Ca1 + Ca2 S_bar + Ca3 W_bar)
 *   S_bar = (k/eps) sqrt(S2)
 *   W_bar = (k/eps) sqrt(W2)
 */
Real Cmu_nonlinear(const Real Ca0,
                   const Real Ca1,
                   const Real Ca2,
                   const Real Ca3,
                   const StrainRotationInvariants & inv,
                   const Real k,
                   const Real eps);

/**
 * Turbulent production G_k with optional compressibility terms.
 *
 * Gk = mu_t S^2
 *    - 2/3 rho k div(u) (if include_compressibility_terms)
 *    - 2/3 mu_t (div(u))^2 (if include_compressibility_terms)
 */
Real computeGk(const Real mu_t,
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
Real computeGb(const Real beta,
               const Real mu_t,
               const Real Pr_t,
               const libMesh::VectorValue<Real> & grad_T,
               const libMesh::VectorValue<Real> & g);

/**
 * Compressibility modification gamma_M.
 *
 * gamma_M = rho C_M k epsilon / c^2
 */
Real computeGammaM(const Real rho, const Real C_M, const Real k, const Real eps, const Real c);

/**
 * Yap correction gamma_y.
 *
 * gamma_y = C_w epsilon^2 / k * max[(l / l_eps - 1) (l / l_eps)^2, 0]
 */
Real computeGammaY(const Real C_w, const Real eps, const Real k, const Real l, const Real l_eps);

/**
 * Low-Re additional production G' for the standard k-epsilon model.
 *
 * G' = D f2 (Gk + 2 mu_t k / d^2) exp(-E Re_d^2)
 */
Real computeGprime(const Real D,
                   const Real E,
                   const Real f2,
                   const Real Gk,
                   const Real mu_t,
                   const Real k,
                   const Real d,
                   const Real Re_d);

/**
 * Utility function to compute velocity gradient tensor in multiple dimensions.
 *
 * This function is dimension-agnostic:
 *  - v == nullptr implies 1D,
 *  - w == nullptr and v != nullptr implies 2D,
 *  - v != nullptr and w != nullptr implies 3D.
 */
RankTwoTensor computeVelocityGradient(const Moose::Functor<Real> & u,
                                      const Moose::Functor<Real> * v,
                                      const Moose::Functor<Real> * w,
                                      const Moose::ElemArg & elem_arg,
                                      const Moose::StateArg & state);

/**
 * Build the non-linear Reynolds stress contribution T_RANS,NL for the
 * Standard k-epsilon model using either the quadratic or cubic
 * constitutive relation.
 *
 * S and W are built from grad_u:
 *   S = 0.5 (grad_u + grad_u^T)
 *   W = 0.5 (grad_u - grad_u^T)
 */
RankTwoTensor computeTRANS_NL(const NonlinearConstitutiveRelation model,
                              const RankTwoTensor & grad_u,
                              const StrainRotationInvariants & inv,
                              const Real mu_t,
                              const Real k,
                              const Real eps);

/**
 * Convenience helper to turn T_RANS,NL into the scalar non-linear
 * production term
 *
 *   G_nl = T_RANS,NL : grad(u)
 */
Real computeGnl(const NonlinearConstitutiveRelation model,
                const RankTwoTensor & grad_u,
                const StrainRotationInvariants & inv,
                const Real mu_t,
                const Real k,
                const Real eps);

/**
 * Curvature / rotation correction factor f_c, based on the local
 * strain- and rotation-rate invariants.
 *
 * This is used to scale the production in the Realizable k-epsilon
 * variants:
 *   P_k = f_c G_k
 *   P_eps = f_c S_k
 */
Real computeCurvatureFactor(const CurvatureCorrectionModel model,
                            const StrainRotationInvariants & inv);

/**
 * Gradient method selector for turbulence production terms.
 *
 * MooseFunctor  – use the gradient cached by MOOSE (Green-Gauss by default).
 * LocalLeastSquares – reconstruct a local inverse-distance-weighted least-squares
 *                     gradient from face-adjacent cell values.  Bypasses the MOOSE
 *                     gradient cache entirely and is therefore insensitive to the
 *                     gradient method configured at the variable level.  Falls back
 *                     to the MOOSE functor gradient when the local stencil is
 *                     degenerate (boundary cell with too few neighbors, singular LS
 *                     matrix, etc.).
 */
enum class TurbVelocityGradientMethod
{
  MooseFunctor,
  LocalLeastSquares
};

/**
 * Local inverse-distance-weighted least-squares velocity gradient.
 *
 * For each face-adjacent neighbor N of element P:
 *   (u_N - u_P) ≈ ∇u · (x_N - x_P)
 *
 * The overdetermined system is solved in the least-squares sense with
 * weights w = 1/|x_N - x_P|^2.  The analytical 2×2 / 3×3 matrix
 * inverse is used so that no external dense-algebra library is needed.
 *
 * Falls back to the MOOSE functor gradient when:
 *   - the element has no interior (non-boundary) neighbors, or
 *   - the resulting normal-equation matrix is nearly singular.
 *
 * @param u_var   x-velocity functor
 * @param v_var   y-velocity functor (nullptr in 1D)
 * @param w_var   z-velocity functor (nullptr in 1D/2D)
 * @param elem    raw libMesh element (needed for neighbor traversal)
 * @param state   solution state
 */
RankTwoTensor computeVelocityGradientLS(const Moose::Functor<Real> & u_var,
                                        const Moose::Functor<Real> * v_var,
                                        const Moose::Functor<Real> * w_var,
                                        const Elem * elem,
                                        const Moose::StateArg & state);

/**
 * Compute strain/rotation invariants using the chosen gradient method.
 *
 * When method == LocalLeastSquares the velocity gradient is computed by
 * computeVelocityGradientLS; otherwise computeVelocityGradient is used.
 */
StrainRotationInvariants
computeStrainRotationInvariantsEx(const Moose::Functor<Real> & u,
                                  const Moose::Functor<Real> * v,
                                  const Moose::Functor<Real> * w,
                                  const Elem * elem,
                                  const Moose::ElemArg & elem_arg,
                                  const Moose::StateArg & state,
                                  TurbVelocityGradientMethod method);

/**
 * Kato–Launder (1993) modified shear production.
 *
 *   G_k^{KL} = μ_t |S| |Ω|
 *            = μ_t √S2 · √W2
 *
 * This eliminates the stagnation-point anomaly: in regions of pure
 * irrotational strain |Ω| → 0, so production → 0, preventing the
 * unphysical k blow-up at stagnation points.
 *
 * Optional compressibility corrections are applied identically to computeGk.
 */
Real computeGkKatoLaunder(const Real mu_t,
                           const Real S2,
                           const Real W2,
                           const Real rho,
                           const Real k,
                           const Real div_u,
                           const bool include_compressibility_terms);

/**
 * Turbulent time scale with Kolmogorov lower bound (Durbin 1996).
 *
 *   T = max(k/ε,  C_t √(ν/ε))
 *
 * The Kolmogorov lower bound prevents the time scale from collapsing to
 * zero as ε → ∞ in near-wall / highly dissipative cells, which would
 * otherwise make the implicit destruction coefficient in the ε-equation
 * arbitrarily large and cause convergence failure.
 *
 * @param k    turbulent kinetic energy
 * @param eps  turbulent dissipation rate
 * @param nu   kinematic viscosity
 * @param Ct   model constant (default 6.0 in STAR-CCM+)
 */
Real limitTurbTimeScale(const Real k, const Real eps, const Real nu, const Real Ct);

} // namespace NS
