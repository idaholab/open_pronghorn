//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TurbulenceMethods.h"

namespace NS
{

StrainRotationInvariants
computeStrainRotationInvariants(const Moose::Functor<Real> & u,
                                const Moose::Functor<Real> * v,
                                const Moose::Functor<Real> * w,
                                const Moose::ElemArg & elem_arg,
                                const Moose::StateArg & state)
{
  const auto & grad_u = u.gradient(elem_arg, state);

  // du/dx components
  const Real dux = grad_u(0);
  Real duy = 0.0;
  Real duz = 0.0;

  // dv/dx components
  Real dvx = 0.0, dvy = 0.0, dvz = 0.0;
  // dw/dx components
  Real dwx = 0.0, dwy = 0.0, dwz = 0.0;

  if (v)
  {
    const auto & grad_v = v->gradient(elem_arg, state);
    dvx = grad_v(0);
    dvy = grad_v(1);
    if (w)
      dvz = grad_v(2);
  }

  if (w)
  {
    const auto & grad_w = w->gradient(elem_arg, state);
    dwx = grad_w(0);
    dwy = grad_w(1);
    dwz = grad_w(2);
  }

  if (v)
  {
    duy = grad_u(1);
    if (w)
      duz = grad_u(2);
  }

  // Divergence
  Real div_u = dux;
  if (v)
    div_u += dvy;
  if (w)
    div_u += dwz;

  // Symmetric strain-rate tensor S_ij
  const Real Sxx = dux;
  const Real Syy = v ? dvy : 0.0;
  const Real Szz = w ? dwz : 0.0;
  const Real Sxy = v ? 0.5 * (duy + dvx) : 0.0;
  const Real Sxz = w ? 0.5 * (duz + dwx) : 0.0;
  const Real Syz = (v && w) ? 0.5 * (dvz + dwy) : 0.0;

  // Anti-symmetric rotation tensor W_ij
  const Real Wxy = v ? 0.5 * (duy - dvx) : 0.0;
  const Real Wxz = w ? 0.5 * (duz - dwx) : 0.0;
  const Real Wyz = (v && w) ? 0.5 * (dvz - dwy) : 0.0;

  // Invariants (2 S_ij S_ij and 2 W_ij W_ij)
  const Real S2 = 2.0 *
                  (Sxx * Sxx + Syy * Syy + Szz * Szz +
                   2.0 * (Sxy * Sxy + Sxz * Sxz + Syz * Syz));

  const Real W2 = 4.0 * (Wxy * Wxy + Wxz * Wxz + Wyz * Wyz);

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
  out.mu_ratio = 0.42 * std::pow(Re_d, 0.25) *
                 (1.0 - std::exp(-Re_d / 70.0));
  return out;
}

TwoLayerLengths
twoLayerNorrisReynolds(const Real Cmu, const Real d, const Real Re_d)
{
  const Real cl = cl_from_Cmu(Cmu);
  TwoLayerLengths out;
  out.l_eps = cl * d * Re_d / (Re_d + 5.3);
  out.mu_ratio = 0.42 * std::pow(Re_d, 0.25) *
                 (1.0 - std::exp(-Re_d / 50.5));
  return out;
}

TwoLayerLengths
twoLayerXu(const Real /*Cmu*/, const Real d, const Real /*Re_d*/, const Real yv_star)
{
  TwoLayerLengths out;
  out.l_eps = 8.8 * d /
              (1.0 + 10.0 / yv_star + 5.15e-2 * yv_star);
  out.mu_ratio = 0.544 * yv_star /
                 (1.0 + 5.025e-4 * std::pow(yv_star, 1.65));
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
  // Gb = beta mu_t / Pr_t * (gradT · g)
  return beta * mu_t / Pr_t * (grad_T * g);
}

Real
computeGammaM(const Real rho,
              const Real C_M,
              const Real k,
              const Real eps,
              const Real c)
{
  // gamma_M = rho C_M k epsilon / c^2
  return rho * C_M * k * eps / (c * c);
}

Real
computeGammaY(const Real C_w,
              const Real eps,
              const Real k,
              const Real l,
              const Real l_eps)
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

} // namespace NS
