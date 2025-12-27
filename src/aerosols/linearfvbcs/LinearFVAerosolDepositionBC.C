//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVAerosolDepositionBC.h"

registerMooseObject("OpenPronghornApp", LinearFVAerosolDepositionBC);

InputParameters
LinearFVAerosolDepositionBC::validParams()
{
  // Start from the generic Robin FV advection-diffusion BC
  InputParameters params = LinearFVAdvectionDiffusionFunctorRobinBC::validParams();

  params.addClassDescription(
      "FV Robin boundary condition for aerosol deposition at walls. "
      "Computes Lai & Nazaroff (2000) deposition velocity internally "
      "using friction velocity, Brownian diffusion, and gravitational settling. "
      "Imposes grad(C).n + (V_d/D) C = 0 so that the diffusive flux is V_d * C.");

  // Fluid properties
  params.addRequiredParam<MooseFunctorName>(
      "rho", "Carrier fluid density functor [kg/m^3] for the aerosol carrier gas.");
  params.addRequiredParam<MooseFunctorName>(
      "mu", "Carrier fluid dynamic viscosity functor [Pa.s] for the carrier gas.");
  params.addRequiredParam<MooseFunctorName>("T", "Carrier fluid temperature functor [K].");

  // Friction velocity at the wall
  params.addRequiredParam<MooseFunctorName>("u_star",
                                            "Friction velocity functor u_* [m/s] at the wall.");

  // Particle properties
  params.addRequiredParam<Real>("particle_diameter", "Aerosol particle diameter [m].");
  params.addRequiredParam<Real>("particle_density", "Aerosol particle density [kg/m^3].");

  // Gas mean free path (for Cunningham slip and Brownian diffusivity)
  params.addParam<Real>("mean_free_path",
                        0.066e-6, // 0.066 micron, typical air value
                        "Mean free path of the carrier gas [m].");

  // Gravity
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0.0, 0.0, -9.81), "Gravitational acceleration vector [m/s^2].");

  // Diffusion coefficient (must match the coefficient used in the PDE)
  params.addParam<MooseFunctorName>("diffusion_coeff", 1.0, "The diffusion coefficient.");

  return params;
}

LinearFVAerosolDepositionBC::LinearFVAerosolDepositionBC(const InputParameters & parameters)
  : LinearFVAdvectionDiffusionFunctorRobinBC(parameters),
    _rho(getFunctor<Real>("rho")),
    _mu(getFunctor<Real>("mu")),
    _T(getFunctor<Real>("T")),
    _u_star(getFunctor<Real>("u_star")),
    _d_p(getParam<Real>("particle_diameter")),
    _rho_p(getParam<Real>("particle_density")),
    _lambda_mfp(getParam<Real>("mean_free_path")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _diffusion_coeff(getFunctor<Real>("diffusion_coeff"))
{
}

// -----------------------------------------------------------------------------
// Robin coefficients
// -----------------------------------------------------------------------------

Real
LinearFVAerosolDepositionBC::getAlpha(Moose::FaceArg /*face*/, Moose::StateArg /*state*/) const
{
  // We enforce: grad(C).n + (V_d / D) C = 0
  // => alpha = 1, beta = V_d / D, gamma = 0
  return 1.0;
}

Real
LinearFVAerosolDepositionBC::getBeta(Moose::FaceArg /*face*/, Moose::StateArg /*state*/) const
{
  const auto state = determineState();

  // Use collocated face arguments for the diffusion coefficient and deposition velocity,
  // consistent with the previous implementation.
  const auto diff_face = makeCDFace(*_current_face_info);
  const Real D = _diffusion_coeff(diff_face, state);

  if (D <= 0.0)
    return 0.0;

  const auto dep_face = singleSidedFaceArg(_current_face_info);
  const Real Vd = computeDepositionVelocity(dep_face);

  if (Vd <= 0.0)
    return 0.0;

  // beta = V_d / D so that -D * grad(C).n = V_d * C
  return Vd / D;
}

Real
LinearFVAerosolDepositionBC::getGamma(Moose::FaceArg /*face*/, Moose::StateArg /*state*/) const
{
  // No inhomogeneous term: flux is purely proportional to C
  return 0.0;
}

// -----------------------------------------------------------------------------
// Deposition-velocity model (unchanged Lai & Nazaroff implementation)
// -----------------------------------------------------------------------------

Real
LinearFVAerosolDepositionBC::computeDepositionVelocity(const Moose::FaceArg & face) const
{
  const auto state = determineState();

  const Real mu_face = _mu(face, state);
  const Real rho_face = _rho(face, state);
  const Real T_face = _T(face, state);
  const Real u_star = _u_star(face, state);

  if (mu_face <= 0.0 || rho_face <= 0.0 || T_face <= 0.0 || u_star <= 0.0)
    return 0.0;

  const auto & normal = _current_face_info->normal();
  const Real gmag = _gravity.norm();

  // Gravitational settling velocity vector at the face
  const RealVectorValue v_gs_vec = computeGravitationalSettlingVelocity(mu_face, rho_face);

  // Wall-normal component magnitude of settling velocity
  const Real Vgs = std::fabs(v_gs_vec * normal);

  // Classify surface orientation using n . g
  Real Vd = 0.0;

  if (gmag > 0.0)
  {
    const Real alignment = (normal * _gravity) / gmag;

    // Mostly horizontal: |alignment| ~ 1 -> floor or ceiling
    const Real horiz_cutoff = 0.5; // cos(60°), tolerate some slope

    if (std::fabs(alignment) > horiz_cutoff)
    {
      // normal // gravity -> floor; normal anti-parallel to gravity -> ceiling
      if (alignment > 0.0)
        Vd = computeFloorDepositionVelocity(Vgs, u_star);
      else
        Vd = computeCeilingDepositionVelocity(Vgs, u_star);
    }
    else
    {
      // Vertical surface
      Vd = computeVerticalDepositionVelocity(u_star, mu_face, rho_face, T_face);
    }
  }
  else
  {
    // No gravity direction (unusual); fall back to vertical-wall model
    Vd = computeVerticalDepositionVelocity(u_star, mu_face, rho_face, T_face);
  }

  if (!std::isfinite(Vd) || Vd < 0.0)
    Vd = 0.0;

  return Vd;
}

Real
LinearFVAerosolDepositionBC::computeCeilingDepositionVelocity(Real Vgs, Real u_star) const
{
  if (u_star <= 0.0)
    return 0.0;

  const Real ratio = Vgs / u_star;

  // Small-argument limit: V_d → u_*
  if (std::fabs(ratio) < 1e-8)
    return u_star;

  const Real denom = 1.0 - std::exp(-ratio);
  if (std::fabs(denom) < 1e-12)
    return u_star;

  return Vgs / denom;
}

Real
LinearFVAerosolDepositionBC::computeFloorDepositionVelocity(Real Vgs, Real u_star) const
{
  if (u_star <= 0.0)
    return 0.0;

  const Real ratio = Vgs / u_star;

  // Small-argument limit: V_d → u_*
  if (std::fabs(ratio) < 1e-8)
    return u_star;

  const Real denom = std::exp(ratio) - 1.0;
  if (std::fabs(denom) < 1e-12)
    return u_star;

  return Vgs / denom;
}

Real
LinearFVAerosolDepositionBC::computeVerticalDepositionVelocity(Real u_star,
                                                               Real mu,
                                                               Real rho,
                                                               Real T) const
{
  if (u_star <= 0.0)
    return 0.0;

  const Real I = computeIntegralI(u_star, mu, rho, T);
  if (I <= 0.0)
    return 0.0;

  return u_star / I;
}

Real
LinearFVAerosolDepositionBC::computeIntegralI(Real u_star, Real mu, Real rho, Real T) const
{
  if (mu <= 0.0 || rho <= 0.0 || T <= 0.0 || u_star <= 0.0 || _d_p <= 0.0)
    return 0.0;

  // Brownian diffusivity (slip-corrected) for aerosol
  const Real D_B = computeBrownianDiffusivity(mu, T);
  if (D_B <= 0.0)
    return 0.0;

  // Kinematic viscosity
  const Real nu = mu / rho;

  // Schmidt number Sc = nu / D_B
  const Real Sc = nu / D_B;
  if (Sc <= 0.0)
    return 0.0;

  const Real Sc_inv = 1.0 / Sc;
  const Real Sc_m13 = std::pow(Sc, -1.0 / 3.0);

  // Dimensionless wall coordinate r+ = d_p u_* / (2 nu)
  const Real r_plus = _d_p * u_star / (2.0 * nu);

  // Lai & Nazaroff auxiliary functions a and b (Table S5)
  const Real A0 = 10.92 * Sc_m13;

  // a: evaluated at y+ = 4.3
  const Real y0_plus = 4.3;
  const Real num_a = std::pow(A0 + y0_plus, 3.0);
  const Real den_a = Sc_inv + 0.0609;
  const Real log_a = 0.5 * std::log(num_a / den_a);
  const Real atan_a = std::sqrt(3.0) * std::atan((2.0 * y0_plus - A0) / (std::sqrt(3.0) * A0));
  const Real a = log_a + atan_a;

  // b: evaluated at r+
  const Real num_b = std::pow(A0 + r_plus, 3.0);
  const Real den_b = Sc_inv + 7.669e-4 * std::pow(r_plus, 3.0);
  const Real log_b = 0.5 * std::log(num_b / den_b);
  const Real atan_b = std::sqrt(3.0) * std::atan((2.0 * r_plus - A0) / (std::sqrt(3.0) * A0));
  const Real b = log_b + atan_b;

  // Integral I
  Real I = 3.64 * std::pow(Sc, 2.0 / 3.0) * (a - b) + 39.0;

  if (!std::isfinite(I) || I <= 0.0)
    I = 0.0;

  return I;
}

Real
LinearFVAerosolDepositionBC::computeCunninghamSlipCorrection() const
{
  if (_d_p <= 0.0 || _lambda_mfp <= 0.0)
    return 1.0;

  const Real Kn_ratio = _lambda_mfp / _d_p;
  return 1.0 + Kn_ratio * (2.34 + 1.05 * std::exp(-0.39 * _d_p / _lambda_mfp));
}

Real
LinearFVAerosolDepositionBC::computeBrownianDiffusivity(Real mu, Real T) const
{
  if (mu <= 0.0 || T <= 0.0 || _d_p <= 0.0)
    return 0.0;

  // Boltzmann constant [J/K]
  static const Real k_B = 1.380649e-23;

  const Real C_c = computeCunninghamSlipCorrection();

  // Stokes–Einstein with slip correction: D_B = k_B T C_c / (3 π μ d_p)
  const Real D_B = k_B * T * C_c / (3.0 * libMesh::pi * mu * _d_p);

  return D_B > 0.0 && std::isfinite(D_B) ? D_B : 0.0;
}

RealVectorValue
LinearFVAerosolDepositionBC::computeGravitationalSettlingVelocity(Real mu, Real /*rho*/) const
{
  if (mu <= 0.0 || _d_p <= 0.0)
    return RealVectorValue(0.0, 0.0, 0.0);

  const Real C_c = computeCunninghamSlipCorrection();

  // Particle relaxation time tau_p = rho_p d_p^2 C_c / (18 μ)
  const Real tau_p = (_rho_p * _d_p * _d_p * C_c) / (18.0 * mu);

  // Settling velocity v_gs = tau_p * g
  return tau_p * _gravity;
}
