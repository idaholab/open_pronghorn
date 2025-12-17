//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVAerosolDriftFlux.h"
#include "LinearFVBoundaryCondition.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", LinearFVAerosolDriftFlux);

InputParameters
LinearFVAerosolDriftFlux::validParams()
{
  auto params = LinearFVFluxKernel::validParams();
  params.addClassDescription(
      "Adds the aerosol drift-flux contribution (rho * U_drift . n * C_face) "
      "to a transported scalar in the linear finite-volume discretization. "
      "The drift velocity is computed internally from gravitational settling "
      "and thermophoresis.");

  params += Moose::FV::advectedInterpolationParameter();

  // Fluid properties
  params.addRequiredParam<MooseFunctorName>("rho", "Carrier fluid density functor [kg/m^3].");
  params.addRequiredParam<MooseFunctorName>("mu",
                                            "Carrier fluid dynamic viscosity functor [Pa.s].");
  params.addRequiredParam<MooseFunctorName>("T", "Carrier fluid temperature functor [K].");

  // Particle properties
  params.addRequiredParam<Real>("particle_diameter", "Aerosol particle diameter [m].");
  params.addRequiredParam<Real>("particle_density", "Aerosol particle density [kg/m^3].");

  // Gas mean free path (for Cunningham slip, etc.)
  params.addParam<Real>("mean_free_path",
                        0.066e-6, // 0.066 micron
                        "Mean free path of the carrier gas [m] used in the "
                        "Cunningham slip correction.");

  // Thermophoresis options
  params.addParam<Real>("k_th_const",
                        0.5,
                        "Constant thermophoretic force coefficient k_th "
                        "used when use_constant_kth=true.");
  params.addParam<bool>("use_constant_kth",
                        true,
                        "If true, use k_th_const; otherwise, use a simple "
                        "diameter-based approximation for k_th.");

  // Gravity
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0.0, 0.0, -9.81), "Gravitational acceleration vector [m/s^2].");

  params.addParam<bool>("force_boundary_execution",
                        true,
                        "This flux kernel should execute on boundaries by default.");

  MooseEnum coeff_interp_method("average harmonic", "harmonic");
  params.addParam<MooseEnum>("density_interp_method",
                             coeff_interp_method,
                             "Interpolation method for the density at faces.");
  params.addParam<MooseEnum>("viscosity_interp_method",
                             coeff_interp_method,
                             "Interpolation method for the viscosity at faces.");

  return params;
}

LinearFVAerosolDriftFlux::LinearFVAerosolDriftFlux(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _rho(getFunctor<Real>("rho")),
    _mu(getFunctor<Real>("mu")),
    _T(getFunctor<Real>("T")),
    _density_interp_method(
        Moose::FV::selectInterpolationMethod(getParam<MooseEnum>("density_interp_method"))),
    _viscosity_interp_method(
        Moose::FV::selectInterpolationMethod(getParam<MooseEnum>("viscosity_interp_method"))),
    _d_p(getParam<Real>("particle_diameter")),
    _rho_p(getParam<Real>("particle_density")),
    _lambda_mfp(getParam<Real>("mean_free_path")),
    _k_th_const(getParam<Real>("k_th_const")),
    _use_constant_kth(getParam<bool>("use_constant_kth")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _face_mass_flux(0.0),
    _interp_coeffs(std::make_pair(0.0, 0.0))
{
  if (dynamic_cast<const MooseLinearVariableFV<Real> *>(&_T))
    dynamic_cast<MooseLinearVariableFV<Real> *>(&_subproblem.getStandardVariable(_tid, "T"))
        ->computeCellGradients();

  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");
}

void
LinearFVAerosolDriftFlux::computeFlux()
{
  const auto & normal = _current_face_info->normal();
  const auto state = determineState();
  const bool on_boundary = Moose::FV::onBoundary(*this, *_current_face_info);

  Moose::FaceArg face_arg;
  if (on_boundary)
    face_arg = singleSidedFaceArg(_current_face_info);
  else
    face_arg = makeCDFace(*_current_face_info);

  // Interpolate density and viscosity at the face
  Real rho_face;
  Real mu_face;

  if (on_boundary)
  {
    rho_face = _rho(face_arg, state);
    mu_face = _mu(face_arg, state);
  }
  else
  {
    const auto elem_arg = makeElemArg(_current_face_info->elemPtr());
    const auto neigh_arg = makeElemArg(_current_face_info->neighborPtr());

    Moose::FV::interpolate(_density_interp_method,
                           rho_face,
                           _rho(elem_arg, state),
                           _rho(neigh_arg, state),
                           *_current_face_info,
                           true);

    Moose::FV::interpolate(_viscosity_interp_method,
                           mu_face,
                           _mu(elem_arg, state),
                           _mu(neigh_arg, state),
                           *_current_face_info,
                           true);
  }

  // Compute individual drift components at the face
  const RealVectorValue v_g = computeGravitationalDrift(face_arg, mu_face, rho_face);
  const RealVectorValue v_th = computeThermophoreticDrift(face_arg, mu_face, rho_face);

  // Total drift velocity at the face
  const RealVectorValue U_drift = v_g + v_th;

  const Real udotn = normal * U_drift;

  // Drift mass flux through the face: m_drift = rho_face * (U_drift.n)
  _face_mass_flux = rho_face * udotn;
}

RealVectorValue
LinearFVAerosolDriftFlux::computeGravitationalDrift(const Moose::FaceArg & /*face_arg*/,
                                                    Real mu_face,
                                                    Real /*rho_face*/) const
{
  const auto state = determineState();

  // Cunningham slip correction factor
  const Real Kn_ratio = _d_p > 0.0 ? _lambda_mfp / _d_p : 0.0;
  const Real C_mu = 1.0 + Kn_ratio * (2.34 + 1.05 * std::exp(-0.39 * _d_p / _lambda_mfp));

  // Particle relaxation time: tau_p = rho_p * d_p^2 * C_mu / (18 * mu)
  const Real tau_p = (_rho_p * _d_p * _d_p * C_mu) / (18.0 * mu_face);

  // Settling velocity vector: v_gs = tau_p * g
  return tau_p * _gravity;
}

RealVectorValue
LinearFVAerosolDriftFlux::computeThermophoreticDrift(const Moose::FaceArg & face_arg,
                                                     Real mu_face,
                                                     Real rho_face) const
{
  const auto state = determineState();

  const Real T_face = _T(face_arg, state);
  if (T_face == 0.0)
    return RealVectorValue(0.0, 0.0, 0.0);

  const auto grad_T_raw = _T.gradient(face_arg, state);
  libMesh::VectorValue<Real> gradT_face(grad_T_raw(0), 0., 0.);
  if (_dim >= 2)
    gradT_face(1) = grad_T_raw(1);
  if (_dim >= 3)
    gradT_face(2) = grad_T_raw(2);

  // Thermophoretic coefficient k_th:
  Real k_th = _k_th_const;
  if (!_use_constant_kth)
  {
    // Simple diameter-based approximation:
    // for small particles (d_p < lambda) we take k_th ≈ 0.5; for larger particles
    // you can refine this correlation as needed.
    if (_d_p <= _lambda_mfp)
      k_th = 0.5;
    else
      k_th = 0.5; // placeholder for a more detailed correlation
  }

  // Thermophoretic drift velocity:
  // v_th = - k_th * nu / (rho * T) * gradT
  const Real factor = -k_th * mu_face / (rho_face * T_face);

  return factor * gradT_face;
}

Real
LinearFVAerosolDriftFlux::computeElemMatrixContribution()
{
  // Linear advective-like drift flux contribution on the element side:
  //   F_face = m_drift * C_face
  //   => dF/dC_elem = interp_coeffs.first * m_drift
  return _interp_coeffs.first * _face_mass_flux * _current_face_area;
}

Real
LinearFVAerosolDriftFlux::computeNeighborMatrixContribution()
{
  // Same as above, but for the neighbor side
  return _interp_coeffs.second * _face_mass_flux * _current_face_area;
}

Real
LinearFVAerosolDriftFlux::computeElemRightHandSideContribution()
{
  // Flux is linear in C, so no explicit RHS term is needed here.
  return 0.0;
}

Real
LinearFVAerosolDriftFlux::computeNeighborRightHandSideContribution()
{
  // Flux is linear in C, so no explicit RHS term is needed here.
  return 0.0;
}

Real
LinearFVAerosolDriftFlux::computeBoundaryRHSContribution(const LinearFVBoundaryCondition & /*bc*/)
{
  // Currently treat boundary faces the same way as internal faces via
  // matrix contributions; no extra explicit RHS term.
  return 0.0;
}

void
LinearFVAerosolDriftFlux::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  // Compute drift mass flux for the current face
  computeFlux();

  // Compute interpolation coefficients for the scalar at the face using the
  // sign of the drift mass flux for upwinding/averaging.
  _interp_coeffs =
      interpCoeffs(_advected_interp_method, *_current_face_info, true, _face_mass_flux);
}
