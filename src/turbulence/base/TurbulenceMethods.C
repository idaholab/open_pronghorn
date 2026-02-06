//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TurbulenceMethods.h"
#include "SystemBase.h"

namespace NS
{

StrainRotationInvariants
computeStrainRotationInvariants(const Moose::Functor<Real> & u,
                                const Moose::Functor<Real> * v,
                                const Moose::Functor<Real> * w,
                                const Moose::ElemArg & elem_arg,
                                const Moose::StateArg & state)
{
  // Reuse the common gradient builder
  const RankTwoTensor grad_vel = computeVelocityGradient(u, v, w, elem_arg, state);

  // Divergence = tr(grad u)
  const Real div_u = grad_vel(0, 0) + grad_vel(1, 1) + grad_vel(2, 2);

  // Symmetric strain-rate tensor and anti-symmetric rotation tensor
  const RankTwoTensor S = 0.5 * (grad_vel + grad_vel.transpose());
  const RankTwoTensor W = 0.5 * (grad_vel - grad_vel.transpose());

  // Invariants: S2 = 2 S_ij S_ij, W2 = 2 W_ij W_ij
  const Real S2 = 2.0 * S.doubleContraction(S);
  const Real W2 = 2.0 * W.doubleContraction(W);

  StrainRotationInvariants out;
  out.S2 = S2;
  out.W2 = W2;
  out.div_u = div_u;
  return out;
}

Real
cl_from_Cmu(const Real Cmu)
{
  return 0.42 * std::pow(Cmu, -0.75);
}

TwoLayerLengths
twoLayerWolfstein(const Real Cmu, const Real d, const Real Re_d)
{
  const Real cl = cl_from_Cmu(Cmu);
  TwoLayerLengths out;
  out.l_eps = cl * d * (1.0 - std::exp(-Re_d / (2.0 * cl)));
  out.mu_ratio = 0.42 * std::pow(Re_d, 0.25) * (1.0 - std::exp(-Re_d / 70.0));
  return out;
}

TwoLayerLengths
twoLayerNorrisReynolds(const Real Cmu, const Real d, const Real Re_d)
{
  const Real cl = cl_from_Cmu(Cmu);
  TwoLayerLengths out;
  out.l_eps = cl * d * Re_d / (Re_d + 5.3);
  out.mu_ratio = 0.42 * std::pow(Re_d, 0.25) * (1.0 - std::exp(-Re_d / 50.5));
  return out;
}

TwoLayerLengths
twoLayerXu(const Real /*Cmu*/, const Real d, const Real /*Re_d*/, const Real yv_star)
{
  TwoLayerLengths out;
  out.l_eps = 8.8 * d / (1.0 + 10.0 / yv_star + 5.15e-2 * yv_star);
  out.mu_ratio = 0.544 * yv_star / (1.0 + 5.025e-4 * std::pow(yv_star, 1.65));
  return out;
}

Real
f2_SKE_LRe(const Real C, const Real Re_t)
{
  return 1.0 - C * std::exp(-Re_t * Re_t);
}

Real
fmu_SKE_LRe(const Real Cd0, const Real Cd1, const Real Cd2, const Real Re_d)
{
  const Real arg = Cd0 * std::sqrt(Re_d) + Cd1 * Re_d + Cd2 * Re_d * Re_d;
  return 1.0 - std::exp(-arg);
}

Real
f2_RKE(const Real k, const Real nu, const Real eps)
{
  return k / (k + std::sqrt(nu * eps));
}

Real
Cmu_realizable(const Real Ca0,
               const Real Ca1,
               const Real Ca2,
               const Real Ca3,
               const Real S2,
               const Real W2,
               const Real k,
               const Real eps)
{
  const Real k_over_eps = k / eps;
  const Real S_bar = k_over_eps * std::sqrt(S2);
  const Real W_bar = k_over_eps * std::sqrt(W2);
  return Ca0 / (Ca1 + Ca2 * S_bar + Ca3 * W_bar);
}

Real
Cmu_nonlinear(const Real Ca0,
              const Real Ca1,
              const Real Ca2,
              const Real Ca3,
              const StrainRotationInvariants & inv,
              const Real k,
              const Real eps)
{
  if (eps <= 0.0 || k <= 0.0)
    return 0.0;

  const Real k_over_eps = k / eps;
  const Real S_bar = k_over_eps * std::sqrt(std::max(inv.S2, 0.0));
  const Real W_bar = k_over_eps * std::sqrt(std::max(inv.W2, 0.0));

  return Ca0 / (Ca1 + Ca2 * S_bar + Ca3 * W_bar);
}

Real
computeGk(const Real mu_t,
          const Real S2,
          const Real rho,
          const Real k,
          const Real div_u,
          const bool include_compressibility_terms)
{
  Real Gk = mu_t * S2;
  if (include_compressibility_terms)
    Gk -= (2.0 / 3.0) * (rho * k * div_u + mu_t * div_u * div_u);
  return Gk;
}

Real
computeGb(const Real beta,
          const Real mu_t,
          const Real Pr_t,
          const libMesh::VectorValue<Real> & grad_T,
          const libMesh::VectorValue<Real> & g)
{
  return beta * mu_t / Pr_t * (grad_T * g);
}

Real
computeGammaM(const Real rho, const Real C_M, const Real k, const Real eps, const Real c)
{
  return rho * C_M * k * eps / (c * c);
}

Real
computeGammaY(const Real C_w, const Real eps, const Real k, const Real l, const Real l_eps)
{
  if (k <= 0.0)
    return 0.0;

  const Real ratio = l / l_eps;
  const Real tmp = (ratio - 1.0) * ratio * ratio;
  const Real val = std::max(tmp, 0.0);
  return C_w * eps * eps / k * val;
}

Real
computeGprime(const Real D,
              const Real E,
              const Real f2,
              const Real Gk,
              const Real mu_t,
              const Real k,
              const Real d,
              const Real Re_d)
{
  const Real term = Gk + 2.0 * mu_t * k / (d * d);
  return D * f2 * term * std::exp(-E * Re_d * Re_d);
}

RankTwoTensor
computeVelocityGradient(const Moose::Functor<Real> & u_var,
                        const Moose::Functor<Real> * v_var,
                        const Moose::Functor<Real> * w_var,
                        const Moose::ElemArg & elem_arg,
                        const Moose::StateArg & state)
{

  RankTwoTensor grad_vel; // zero initialization

  // grad_u_ij = du_i / dx_j

  // Dim > 0
  const auto & u_grad = u_var.gradient(elem_arg, state);
  grad_vel(0, 0) = u_grad(0);

  if (v_var) // Dim > 1
  {
    const auto & v_grad = v_var->gradient(elem_arg, state);
    grad_vel(0, 1) = u_grad(1);
    grad_vel(1, 0) = v_grad(0);
    grad_vel(1, 1) = v_grad(1);

    if (w_var) // Dim > 2
    {
      const auto & w_grad = w_var->gradient(elem_arg, state);
      grad_vel(0, 2) = u_grad(2);
      grad_vel(1, 2) = v_grad(2);
      grad_vel(2, 0) = w_grad(0);
      grad_vel(2, 1) = w_grad(1);
      grad_vel(2, 2) = w_grad(2);
    }
  }

  return grad_vel;
}

namespace
{
// Model coefficients for the non-linear constitutive relation.
// Quadratic / cubic model constants:
static constexpr Real CNL1 = 0.75;
static constexpr Real CNL2 = 3.75;
static constexpr Real CNL3 = 4.75;
static constexpr Real CNL4 = -10.0;
static constexpr Real CNL5 = -2.0;
static constexpr Real CNL6 = 1000.0;
static constexpr Real CNL7 = 1.0;

// Coefficients for C_mu expression (same as your doc for Ca0..3)
static constexpr Real Ca0 = 0.667;
static constexpr Real Ca1 = 1.25;
static constexpr Real Ca2 = 1.0;
static constexpr Real Ca3 = 0.9;
} // anonymous namespace (don't want users modifying this).

RankTwoTensor
computeTRANS_NL(const NonlinearConstitutiveRelation model,
                const RankTwoTensor & grad_u,
                const StrainRotationInvariants & inv,
                const Real mu_t,
                const Real k,
                const Real eps)
{
  RankTwoTensor Tnl; // default = 0

  if (model == NonlinearConstitutiveRelation::None || mu_t <= 0.0 || k <= 0.0 || eps <= 0.0)
    return Tnl;

  // Build strain- and rotation-rate tensors from grad(u)
  RankTwoTensor S = 0.5 * (grad_u + grad_u.transpose());
  RankTwoTensor W = 0.5 * (grad_u - grad_u.transpose());

  // Invariants given in StrainRotationInvariants:
  // S2 = 2 S_ij S_ij, W2 = 2 W_ij W_ij
  const Real SdotS = 0.5 * inv.S2;
  const Real WdotW = 0.5 * inv.W2;
  const Real Sstar = std::sqrt(std::max(SdotS, 0.0)); // |S| = sqrt(S_ij S_ij)

  // Variable C_mu
  const Real Cmu = Cmu_nonlinear(Ca0, Ca1, Ca2, Ca3, inv, k, eps);

  if (Cmu <= 0.0)
    return Tnl;

  // Denominator used in C1, C2, C3
  const Real denom = (CNL6 + CNL7 * Sstar) * Cmu;
  if (denom <= 0.0)
    return Tnl;

  const Real C1 = CNL1 / denom;
  const Real C2 = CNL2 / denom;
  const Real C3 = CNL3 / denom;

  // Identity tensor
  RankTwoTensor I;
  I.zero();
  I.addIa(1.0);

  // Quadratic bracket
  const RankTwoTensor SS = S * S;
  const RankTwoTensor WW = W * W;
  const RankTwoTensor WS = W * S;
  const RankTwoTensor SWt = S * W.transpose();

  const RankTwoTensor quad_bracket =
      C1 * (SS - (SdotS / 3.0) * I) + C2 * (WS + SWt) + C3 * (WW - (WdotW / 3.0) * I);

  // T_RANS,quad
  const RankTwoTensor Tquad = -4.0 * mu_t * (k / eps) * quad_bracket;

  if (model == NonlinearConstitutiveRelation::Quadratic)
    return Tquad;

  // Cubic additions
  if (model == NonlinearConstitutiveRelation::Cubic)
  {
    const Real C4 = CNL4 * Cmu * Cmu;
    const Real C5 = CNL5 * Cmu * Cmu;
    const Real k2_over_eps2 = (k * k) / (eps * eps);

    RankTwoTensor SSw = (S * S) * W;
    RankTwoTensor WtSS = W.transpose() * (S * S);
    RankTwoTensor term4 = SSw + WtSS;

    RankTwoTensor term5 = (SdotS - WdotW) * (S - W.transpose());

    RankTwoTensor cubic_bracket = C4 * term4 + C5 * term5;

    Tnl = Tquad - 8.0 * mu_t * k2_over_eps2 * cubic_bracket;
  }

  return Tnl;
}

Real
computeGnl(const NonlinearConstitutiveRelation model,
           const RankTwoTensor & grad_u,
           const StrainRotationInvariants & inv,
           const Real mu_t,
           const Real k,
           const Real eps)
{
  if (model == NonlinearConstitutiveRelation::None)
    return 0.0;

  const RankTwoTensor Tnl = computeTRANS_NL(model, grad_u, inv, mu_t, k, eps);

  // G_nl = T_RANS,NL : grad(u)
  return Tnl.doubleContraction(grad_u);
}

Real
computeCurvatureFactor(const CurvatureCorrectionModel model, const StrainRotationInvariants & inv)
{
  // No correction requested
  if (model == CurvatureCorrectionModel::None)
    return 1.0;

  // Spalart–Shur-style curvature correction using S2 and W2.
  // inv.S2 = 2 S_ij S_ij, inv.W2 = 2 W_ij W_ij.
  const Real S2 = std::max(inv.S2, 0.0);
  const Real W2 = std::max(inv.W2, 0.0);

  if (S2 <= 0.0)
    return 1.0;

  const Real SS = 0.5 * S2;
  const Real WW = 0.5 * W2;

  // Dimensionless measures of rotation vs strain
  const Real r_star = WW / std::max(SS, 1e-16);       // rotation/strain ratio
  const Real r_tilde = WW / std::max(SS + WW, 1e-16); // bounded (0,1)

  // Model coefficients -- better not to change
  static constexpr Real Crot1 = 1.0;
  static constexpr Real Crot2 = 2.0;
  static constexpr Real Crot3 = 1.0;
  static constexpr Real Cmax = 1.25;

  // Spalart–Shur rotation/curvature function (simplified form)
  Real f_rot = (1.0 + Crot1) * (2.0 * r_tilde) / (1.0 + r_tilde) *
                   (1.0 - Crot2 * r_star / (1.0 + Crot3 * r_star)) -
               Crot1;

  // Bound and shift to get a positive correction factor
  f_rot = std::max(f_rot, -1.0);
  const Real f_c = std::min(Cmax, 1.0 + f_rot);

  return f_c;
}

} // namespace NS
