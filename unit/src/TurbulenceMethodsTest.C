//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#include "gtest/gtest.h"

#include "TurbulenceMethods.h"

#include <cmath>

// Tolerance for comparisons against values recomputed from the closure definitions. The
// expressions are evaluated in a different order here than in the implementation, so only
// round-off level agreement is expected.
static constexpr Real tol = 1e-12;

TEST(TurbulenceMethods, clFromCmu)
{
  // c_l = 0.42 Cmu^(-3/4)
  EXPECT_NEAR(NS::cl_from_Cmu(0.09), 0.42 / std::pow(0.09, 0.75), tol);

  // The standard C_mu = 0.09 gives the classical c_l ~ 2.556 of the two-layer models
  EXPECT_NEAR(NS::cl_from_Cmu(0.09), 2.556, 1e-3);

  // c_l decreases with C_mu, since the exponent is negative
  EXPECT_GT(NS::cl_from_Cmu(0.05), NS::cl_from_Cmu(0.09));
}

TEST(TurbulenceMethods, twoLayerWolfstein)
{
  const Real Cmu = 0.09;
  const Real d = 0.002;
  const Real Re_d = 40.0;
  const Real cl = NS::cl_from_Cmu(Cmu);

  const auto out = NS::twoLayerWolfstein(Cmu, d, Re_d);

  // l_eps = c_l d (1 - exp(-Re_d / (2 c_l)))
  EXPECT_NEAR(out.l_eps, cl * d * (1.0 - std::exp(-Re_d / (2.0 * cl))), tol);
  // (mu_t/mu)_2L = 0.42 Re_d^(1/4) (1 - exp(-Re_d / 70))
  EXPECT_NEAR(out.mu_ratio, 0.42 * std::pow(Re_d, 0.25) * (1.0 - std::exp(-Re_d / 70.0)), tol);

  // Both length and viscosity ratio vanish at the wall, where d and Re_d go to zero
  const auto at_wall = NS::twoLayerWolfstein(Cmu, 0.0, 0.0);
  EXPECT_NEAR(at_wall.l_eps, 0.0, tol);
  EXPECT_NEAR(at_wall.mu_ratio, 0.0, tol);

  // l_eps is proportional to the wall distance at fixed Re_d
  const auto doubled = NS::twoLayerWolfstein(Cmu, 2.0 * d, Re_d);
  EXPECT_NEAR(doubled.l_eps, 2.0 * out.l_eps, tol);

  // Far from the wall the exponential saturates and l_eps approaches c_l d
  const auto far = NS::twoLayerWolfstein(Cmu, d, 1e4);
  EXPECT_NEAR(far.l_eps, cl * d, 1e-10);
}

TEST(TurbulenceMethods, twoLayerNorrisReynolds)
{
  const Real Cmu = 0.09;
  const Real d = 0.002;
  const Real Re_d = 40.0;
  const Real cl = NS::cl_from_Cmu(Cmu);

  const auto out = NS::twoLayerNorrisReynolds(Cmu, d, Re_d);

  // l_eps = c_l d Re_d / (Re_d + 5.3)
  EXPECT_NEAR(out.l_eps, cl * d * Re_d / (Re_d + 5.3), tol);
  // (mu_t/mu)_2L = 0.42 Re_d^(1/4) (1 - exp(-Re_d / 50.5))
  EXPECT_NEAR(out.mu_ratio, 0.42 * std::pow(Re_d, 0.25) * (1.0 - std::exp(-Re_d / 50.5)), tol);

  // At the wall both vanish
  const auto at_wall = NS::twoLayerNorrisReynolds(Cmu, 0.0, 0.0);
  EXPECT_NEAR(at_wall.l_eps, 0.0, tol);
  EXPECT_NEAR(at_wall.mu_ratio, 0.0, tol);

  // Far from the wall the rational function saturates and l_eps approaches c_l d. The
  // approach is algebraic, Re_d / (Re_d + 5.3), so the residual at Re_d = 1e8 is ~5e-8
  // relative rather than round-off
  const auto far = NS::twoLayerNorrisReynolds(Cmu, d, 1e8);
  EXPECT_NEAR(far.l_eps, cl * d, 1e-9);

  // Norris-Reynolds damps the viscosity less than Wolfstein, since it uses a shorter
  // exponential scale, 50.5 instead of 70
  const auto wolfstein = NS::twoLayerWolfstein(Cmu, d, Re_d);
  EXPECT_GT(out.mu_ratio, wolfstein.mu_ratio);
}

TEST(TurbulenceMethods, twoLayerXu)
{
  const Real d = 0.002;
  const Real yv_star = 30.0;

  // The Xu model ignores C_mu and Re_d, driving both quantities with yv_star only
  const auto out = NS::twoLayerXu(0.09, d, 40.0, yv_star);
  const auto other_cmu = NS::twoLayerXu(0.05, d, 1234.0, yv_star);
  EXPECT_NEAR(out.l_eps, other_cmu.l_eps, tol);
  EXPECT_NEAR(out.mu_ratio, other_cmu.mu_ratio, tol);

  // l_eps = 8.8 d / (1 + 10 / yv* + 5.15e-2 yv*)
  EXPECT_NEAR(out.l_eps, 8.8 * d / (1.0 + 10.0 / yv_star + 5.15e-2 * yv_star), tol);
  // (mu_t/mu)_2L = 0.544 yv* / (1 + 5.025e-4 yv*^1.65)
  EXPECT_NEAR(out.mu_ratio, 0.544 * yv_star / (1.0 + 5.025e-4 * std::pow(yv_star, 1.65)), tol);

  // l_eps remains proportional to the wall distance
  const auto doubled = NS::twoLayerXu(0.09, 2.0 * d, 40.0, yv_star);
  EXPECT_NEAR(doubled.l_eps, 2.0 * out.l_eps, tol);
}

TEST(TurbulenceMethods, f2LowRe)
{
  // f2 = 1 - C exp(-Re_t^2)
  EXPECT_NEAR(NS::f2_SKE_LRe(1.0, 2.0), 1.0 - std::exp(-4.0), tol);

  // At vanishing turbulent Reynolds number the damping is strongest, f2 = 1 - C
  EXPECT_NEAR(NS::f2_SKE_LRe(0.3, 0.0), 0.7, tol);

  // At high turbulent Reynolds number the correction switches off and f2 tends to 1
  EXPECT_NEAR(NS::f2_SKE_LRe(1.0, 50.0), 1.0, tol);

  // C = 0 disables the damping at any Re_t
  EXPECT_NEAR(NS::f2_SKE_LRe(0.0, 0.5), 1.0, tol);
}

TEST(TurbulenceMethods, fmuLowRe)
{
  const Real Cd0 = 0.091, Cd1 = 0.0042, Cd2 = 0.00011;
  const Real Re_d = 12.0;

  // fmu = 1 - exp(-(Cd0 sqrt(Re_d) + Cd1 Re_d + Cd2 Re_d^2))
  const Real arg = Cd0 * std::sqrt(Re_d) + Cd1 * Re_d + Cd2 * Re_d * Re_d;
  EXPECT_NEAR(NS::fmu_SKE_LRe(Cd0, Cd1, Cd2, Re_d), 1.0 - std::exp(-arg), tol);

  // At the wall the viscosity is fully damped
  EXPECT_NEAR(NS::fmu_SKE_LRe(Cd0, Cd1, Cd2, 0.0), 0.0, tol);

  // Far from the wall the damping switches off
  EXPECT_NEAR(NS::fmu_SKE_LRe(Cd0, Cd1, Cd2, 1e5), 1.0, tol);

  // The damping function is a monotone ramp bounded by 0 and 1
  Real previous = 0.0;
  for (const auto i : make_range(1, 40))
  {
    const Real fmu = NS::fmu_SKE_LRe(Cd0, Cd1, Cd2, static_cast<Real>(i));
    EXPECT_GE(fmu, previous);
    EXPECT_GE(fmu, 0.0);
    EXPECT_LE(fmu, 1.0);
    previous = fmu;
  }
}

TEST(TurbulenceMethods, f2Realizable)
{
  const Real k = 0.1, nu = 1e-5, eps = 0.02;

  // f2 = k / (k + sqrt(nu eps))
  EXPECT_NEAR(NS::f2_RKE(k, nu, eps), k / (k + std::sqrt(nu * eps)), tol);

  // The realizable f2 is a fraction below unity, approaching 1 at high turbulence levels
  EXPECT_GT(NS::f2_RKE(k, nu, eps), 0.0);
  EXPECT_LT(NS::f2_RKE(k, nu, eps), 1.0);
  EXPECT_NEAR(NS::f2_RKE(1e8, nu, eps), 1.0, 1e-8);

  // With no viscous contribution the damping disappears exactly
  EXPECT_NEAR(NS::f2_RKE(k, 0.0, eps), 1.0, tol);
}

TEST(TurbulenceMethods, CmuRealizable)
{
  const Real Ca0 = 0.667, Ca1 = 1.25, Ca2 = 1.0, Ca3 = 0.9;
  const Real S2 = 4.0, W2 = 1.0, k = 0.1, eps = 0.05;

  // Cmu = Ca0 / (Ca1 + Ca2 S_bar + Ca3 W_bar), S_bar = (k/eps) sqrt(S2), likewise for W
  const Real S_bar = (k / eps) * std::sqrt(S2);
  const Real W_bar = (k / eps) * std::sqrt(W2);
  EXPECT_NEAR(NS::Cmu_realizable(Ca0, Ca1, Ca2, Ca3, S2, W2, k, eps),
              Ca0 / (Ca1 + Ca2 * S_bar + Ca3 * W_bar),
              tol);

  // Without strain or rotation the expression reduces to Ca0 / Ca1
  EXPECT_NEAR(NS::Cmu_realizable(Ca0, Ca1, Ca2, Ca3, 0.0, 0.0, k, eps), Ca0 / Ca1, tol);

  // Increasing the strain reduces C_mu, which is what bounds the realizable Reynolds stress
  EXPECT_LT(NS::Cmu_realizable(Ca0, Ca1, Ca2, Ca3, 100.0, W2, k, eps),
            NS::Cmu_realizable(Ca0, Ca1, Ca2, Ca3, S2, W2, k, eps));
}

TEST(TurbulenceMethods, CmuNonlinear)
{
  const Real Ca0 = 0.667, Ca1 = 1.25, Ca2 = 1.0, Ca3 = 0.9;
  NS::StrainRotationInvariants inv{/*S2*/ 4.0, /*W2*/ 1.0, /*div_u*/ 0.0};
  const Real k = 0.1, eps = 0.05;

  // Matches the realizable expression when fed the same invariants
  EXPECT_NEAR(NS::Cmu_nonlinear(Ca0, Ca1, Ca2, Ca3, inv, k, eps),
              NS::Cmu_realizable(Ca0, Ca1, Ca2, Ca3, inv.S2, inv.W2, k, eps),
              tol);

  // Degenerate turbulence states return zero rather than dividing by zero
  EXPECT_NEAR(NS::Cmu_nonlinear(Ca0, Ca1, Ca2, Ca3, inv, k, 0.0), 0.0, tol);
  EXPECT_NEAR(NS::Cmu_nonlinear(Ca0, Ca1, Ca2, Ca3, inv, 0.0, eps), 0.0, tol);

  // Negative invariants are clipped to zero before the square roots are taken
  NS::StrainRotationInvariants negative{-1.0, -1.0, 0.0};
  EXPECT_NEAR(NS::Cmu_nonlinear(Ca0, Ca1, Ca2, Ca3, negative, k, eps), Ca0 / Ca1, tol);
}

TEST(TurbulenceMethods, computeGk)
{
  const Real mu_t = 1e-3, S2 = 5.0, rho = 1.2, k = 0.1, div_u = 0.3;

  // Incompressible form: Gk = mu_t S2
  EXPECT_NEAR(NS::computeGk(mu_t, S2, rho, k, div_u, false), mu_t * S2, tol);

  // The divergence is ignored entirely when the compressibility terms are off
  EXPECT_NEAR(NS::computeGk(mu_t, S2, rho, k, 0.0, false),
              NS::computeGk(mu_t, S2, rho, k, div_u, false),
              tol);

  // Compressible form subtracts 2/3 (rho k div_u + mu_t div_u^2)
  EXPECT_NEAR(NS::computeGk(mu_t, S2, rho, k, div_u, true),
              mu_t * S2 - (2.0 / 3.0) * (rho * k * div_u + mu_t * div_u * div_u),
              tol);

  // In an expanding flow the correction reduces the production
  EXPECT_LT(NS::computeGk(mu_t, S2, rho, k, div_u, true),
            NS::computeGk(mu_t, S2, rho, k, div_u, false));

  // At zero divergence both forms agree
  EXPECT_NEAR(
      NS::computeGk(mu_t, S2, rho, k, 0.0, true), NS::computeGk(mu_t, S2, rho, k, 0.0, false), tol);
}

TEST(TurbulenceMethods, computeGb)
{
  const Real beta = 3.4e-3, mu_t = 1e-3, Pr_t = 0.9;
  const libMesh::VectorValue<Real> g(0.0, -9.81, 0.0);

  // Gb = beta mu_t / Pr_t (grad_T . g)
  const libMesh::VectorValue<Real> grad_T(0.0, 5.0, 0.0);
  EXPECT_NEAR(NS::computeGb(beta, mu_t, Pr_t, grad_T, g), beta * mu_t / Pr_t * (grad_T * g), tol);

  // Temperature increasing upwards is a stable stratification: buoyancy destroys turbulence
  EXPECT_LT(NS::computeGb(beta, mu_t, Pr_t, grad_T, g), 0.0);

  // Temperature increasing downwards is unstable: buoyancy produces turbulence
  const libMesh::VectorValue<Real> unstable(0.0, -5.0, 0.0);
  EXPECT_GT(NS::computeGb(beta, mu_t, Pr_t, unstable, g), 0.0);

  // A gradient orthogonal to gravity does not contribute
  const libMesh::VectorValue<Real> horizontal(5.0, 0.0, 0.0);
  EXPECT_NEAR(NS::computeGb(beta, mu_t, Pr_t, horizontal, g), 0.0, tol);
}

TEST(TurbulenceMethods, computeGammaM)
{
  const Real rho = 1.2, C_M = 1.0, k = 0.1, eps = 0.05, c = 340.0;

  // gamma_M = rho C_M k eps / c^2
  EXPECT_NEAR(NS::computeGammaM(rho, C_M, k, eps, c), rho * C_M * k * eps / (c * c), tol);

  // The correction is negligible at low Mach number and grows as the speed of sound drops
  EXPECT_GT(NS::computeGammaM(rho, C_M, k, eps, 30.0), NS::computeGammaM(rho, C_M, k, eps, 340.0));
}

TEST(TurbulenceMethods, computeGammaY)
{
  const Real C_w = 0.83, eps = 0.05, k = 0.1;

  // gamma_y = C_w eps^2 / k max[(l/l_eps - 1)(l/l_eps)^2, 0]
  const Real l = 0.03, l_eps = 0.02;
  const Real ratio = l / l_eps;
  EXPECT_NEAR(NS::computeGammaY(C_w, eps, k, l, l_eps),
              C_w * eps * eps / k * (ratio - 1.0) * ratio * ratio,
              tol);

  // The correction only acts where the length scale exceeds the equilibrium one. Below that
  // the max() clips it to zero rather than letting it change sign
  EXPECT_NEAR(NS::computeGammaY(C_w, eps, k, 0.01, 0.02), 0.0, tol);
  EXPECT_NEAR(NS::computeGammaY(C_w, eps, k, 0.02, 0.02), 0.0, tol);
  EXPECT_GT(NS::computeGammaY(C_w, eps, k, 0.03, 0.02), 0.0);

  // Vanishing turbulent kinetic energy returns zero instead of dividing by zero
  EXPECT_NEAR(NS::computeGammaY(C_w, eps, 0.0, l, l_eps), 0.0, tol);
}

TEST(TurbulenceMethods, computeGprime)
{
  const Real D = 1.0, E = 0.00375, f2 = 0.9, Gk = 2e-3, mu_t = 1e-3, k = 0.1, d = 0.002;
  const Real Re_d = 12.0;

  // G' = D f2 (Gk + 2 mu_t k / d^2) exp(-E Re_d^2)
  EXPECT_NEAR(NS::computeGprime(D, E, f2, Gk, mu_t, k, d, Re_d),
              D * f2 * (Gk + 2.0 * mu_t * k / (d * d)) * std::exp(-E * Re_d * Re_d),
              tol);

  // The extra production is a near-wall term: it decays away from the wall
  EXPECT_LT(NS::computeGprime(D, E, f2, Gk, mu_t, k, d, 100.0),
            NS::computeGprime(D, E, f2, Gk, mu_t, k, d, Re_d));

  // D = 0 disables the term
  EXPECT_NEAR(NS::computeGprime(0.0, E, f2, Gk, mu_t, k, d, Re_d), 0.0, tol);
}

TEST(TurbulenceMethods, computeCurvatureFactor)
{
  NS::StrainRotationInvariants inv{/*S2*/ 4.0, /*W2*/ 1.0, /*div_u*/ 0.0};

  // The 'none' model leaves the production untouched
  EXPECT_NEAR(NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::None, inv), 1.0, tol);

  // Pure strain with no rotation currently returns 0, which suppresses the production
  // entirely. A rotation/curvature correction is expected to leave an irrotational flow
  // untouched, i.e. f_c = 1: with r_star = r_tilde = 0 the Spalart-Shur expression reduces
  // to -Crot1 = -1, and the f_c = 1 + f_rot shift then gives 0. This asserts the behavior
  // as it stands so that a change to the formulation is caught; see the note in the commit
  // message introducing these tests.
  NS::StrainRotationInvariants pure_strain{4.0, 0.0, 0.0};
  EXPECT_NEAR(
      NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::Standard, pure_strain), 0.0, tol);

  // Degenerate state without strain is guarded and falls back to no correction
  NS::StrainRotationInvariants no_strain{0.0, 1.0, 0.0};
  EXPECT_NEAR(
      NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::Standard, no_strain), 1.0, tol);

  // The factor peaks where rotation and strain are comparable, and falls off on either side
  NS::StrainRotationInvariants balanced{4.0, 1.0, 0.0};
  const Real f_c_balanced =
      NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::Standard, balanced);
  EXPECT_GT(f_c_balanced, 0.0);
  EXPECT_GT(f_c_balanced,
            NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::Standard,
                                       NS::StrainRotationInvariants{4.0, 0.1, 0.0}));
  EXPECT_GT(f_c_balanced,
            NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::Standard,
                                       NS::StrainRotationInvariants{4.0, 8.0, 0.0}));

  // The factor stays within the bounds the model imposes, [0, 1.25], over a wide range of
  // rotation-to-strain ratios
  for (const auto i : make_range(1, 60))
  {
    NS::StrainRotationInvariants swept{4.0, 0.2 * static_cast<Real>(i), 0.0};
    const Real f_c = NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::Standard, swept);
    EXPECT_GE(f_c, 0.0);
    EXPECT_LE(f_c, 1.25);
  }

  // Negative invariants are clipped rather than producing a NaN from the square roots
  NS::StrainRotationInvariants negative{-1.0, -1.0, 0.0};
  const Real f_c_negative =
      NS::computeCurvatureFactor(NS::CurvatureCorrectionModel::Standard, negative);
  EXPECT_FALSE(std::isnan(f_c_negative));
  EXPECT_NEAR(f_c_negative, 1.0, tol);
}

TEST(TurbulenceMethods, computeTRANS_NLNoModel)
{
  const RankTwoTensor grad_u(1.0, 0.2, 0.0, 0.3, -1.0, 0.0, 0.0, 0.0, 0.0);
  NS::StrainRotationInvariants inv{4.0, 1.0, 0.0};
  const Real mu_t = 1e-3, k = 0.1, eps = 0.05;

  // The 'none' relation contributes no stress and no production
  const auto none =
      NS::computeTRANS_NL(NS::NonlinearConstitutiveRelation::None, grad_u, inv, mu_t, k, eps);
  EXPECT_NEAR(none.norm(), 0.0, tol);
  EXPECT_NEAR(
      NS::computeGnl(NS::NonlinearConstitutiveRelation::None, grad_u, inv, mu_t, k, eps), 0.0, tol);

  // Degenerate turbulence states are guarded instead of dividing by zero
  for (const auto model :
       {NS::NonlinearConstitutiveRelation::Quadratic, NS::NonlinearConstitutiveRelation::Cubic})
  {
    EXPECT_NEAR(NS::computeTRANS_NL(model, grad_u, inv, 0.0, k, eps).norm(), 0.0, tol);
    EXPECT_NEAR(NS::computeTRANS_NL(model, grad_u, inv, mu_t, 0.0, eps).norm(), 0.0, tol);
    EXPECT_NEAR(NS::computeTRANS_NL(model, grad_u, inv, mu_t, k, 0.0).norm(), 0.0, tol);
  }
}

TEST(TurbulenceMethods, computeTRANS_NLQuadraticAndCubic)
{
  // A pure shear velocity gradient, du/dy = 1
  const RankTwoTensor grad_u(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  // For this gradient S2 = 2 S_ij S_ij = 1 and W2 = 2 W_ij W_ij = 1
  NS::StrainRotationInvariants inv{1.0, 1.0, 0.0};
  const Real mu_t = 1e-3, k = 0.1, eps = 0.05;

  const auto quad =
      NS::computeTRANS_NL(NS::NonlinearConstitutiveRelation::Quadratic, grad_u, inv, mu_t, k, eps);
  const auto cubic =
      NS::computeTRANS_NL(NS::NonlinearConstitutiveRelation::Cubic, grad_u, inv, mu_t, k, eps);

  // Both relations produce a non-zero extra Reynolds stress
  EXPECT_GT(quad.norm(), 0.0);
  EXPECT_GT(cubic.norm(), 0.0);

  // For pure shear both cubic terms vanish identically: S*S*W + W^T*S*S is zero, and
  // S_ij S_ij equals W_ij W_ij so the second term carries a zero prefactor. The cubic
  // relation therefore reduces to the quadratic one in this flow
  EXPECT_NEAR((cubic - quad).norm(), 0.0, 1e-14);

  // The quadratic stress is symmetric for this gradient: it is built from S*S, W*W and the
  // symmetric combination W*S + S*W^T
  EXPECT_NEAR((quad - quad.transpose()).norm(), 0.0, 1e-14);

  // The non-linear stress scales linearly with the turbulent viscosity
  const auto quad_scaled = NS::computeTRANS_NL(
      NS::NonlinearConstitutiveRelation::Quadratic, grad_u, inv, 2.0 * mu_t, k, eps);
  EXPECT_NEAR((quad_scaled - 2.0 * quad).norm(), 0.0, 1e-14);

  // G_nl is the contraction of that stress with the velocity gradient
  EXPECT_NEAR(
      NS::computeGnl(NS::NonlinearConstitutiveRelation::Quadratic, grad_u, inv, mu_t, k, eps),
      quad.doubleContraction(grad_u),
      tol);
  EXPECT_NEAR(NS::computeGnl(NS::NonlinearConstitutiveRelation::Cubic, grad_u, inv, mu_t, k, eps),
              cubic.doubleContraction(grad_u),
              tol);
}

TEST(TurbulenceMethods, computeTRANS_NLGeneralGradient)
{
  // A fully three-dimensional velocity gradient, with the invariants that correspond to it:
  // S2 = 2 S_ij S_ij and W2 = 2 W_ij W_ij
  const RankTwoTensor grad_u(0.5, 1.0, -0.2, 0.3, -0.4, 0.1, 0.7, -0.6, -0.1);
  NS::StrainRotationInvariants inv{/*S2*/ 3.03, /*W2*/ 1.79, /*div_u*/ 0.0};
  const Real mu_t = 1e-3, k = 0.1, eps = 0.05;

  const auto quad =
      NS::computeTRANS_NL(NS::NonlinearConstitutiveRelation::Quadratic, grad_u, inv, mu_t, k, eps);
  const auto cubic =
      NS::computeTRANS_NL(NS::NonlinearConstitutiveRelation::Cubic, grad_u, inv, mu_t, k, eps);

  // Away from pure shear the cubic terms are non-zero, so the two relations differ
  EXPECT_GT((cubic - quad).norm(), 0.0);

  // The quadratic bracket removes the isotropic part of S*S and W*W so the extra stress
  // carries no trace and does not act as a pressure. That holds for the S*S term, but the
  // W*W term is left with a trace: W is antisymmetric, so trace(W*W) = -W_ij W_ij, and
  // subtracting +W_ij W_ij / 3 doubles the trace instead of removing it. The expected value
  // below is 2 * C3 * W_ij W_ij scaled by -4 mu_t k / eps, i.e. the current behavior rather
  // than the trace-free result; see the note in the commit message introducing these tests.
  EXPECT_GT(std::abs(quad.trace()), 1e-6);

  // The S*S contribution on its own is correctly trace free, which isolates the defect to
  // the rotation term: with W = 0 the bracket is trace free as intended
  const RankTwoTensor symmetric_grad(0.5, 0.65, 0.25, 0.65, -0.4, -0.25, 0.25, -0.25, -0.1);
  NS::StrainRotationInvariants sym_inv{/*S2*/ 3.03, /*W2*/ 0.0, /*div_u*/ 0.0};
  const auto sym_quad = NS::computeTRANS_NL(
      NS::NonlinearConstitutiveRelation::Quadratic, symmetric_grad, sym_inv, mu_t, k, eps);
  EXPECT_NEAR(sym_quad.trace(), 0.0, 1e-14);
}
