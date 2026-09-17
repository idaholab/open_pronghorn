#include "MSRLinearFVScalarAdvectionGenericRC.h"

#include "MooseLinearVariableFV.h"
#include "NSFVUtils.h"
#include "NS.h"

registerMooseObject("OpenPronghornApp", MSRLinearFVScalarAdvectionGenericRC);

InputParameters
MSRLinearFVScalarAdvectionGenericRC::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();

  params.addClassDescription(
      "Linear FV passive-scalar advection using the generic RhieChowFaceFluxProvider interface.");

  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "Rhie-Chow face-flux provider used to determine the volumetric face flux.");

  params += Moose::FV::advectedInterpolationParameter();

  params.addParam<MooseFunctorName>("u_slip", "The slip velocity in the x direction.");
  params.addParam<MooseFunctorName>("v_slip", "The slip velocity in the y direction.");
  params.addParam<MooseFunctorName>("w_slip", "The slip velocity in the z direction.");

  return params;
}

MSRLinearFVScalarAdvectionGenericRC::MSRLinearFVScalarAdvectionGenericRC(
    const InputParameters & params)
  : LinearFVFluxKernel(params),
    _mass_flux_provider(
        getUserObject<RhieChowFaceFluxProvider>("rhie_chow_user_object")),
    _advected_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _volumetric_face_flux(0.0),
    _u_slip(isParamValid("u_slip") ? &getFunctor<ADReal>("u_slip") : nullptr),
    _v_slip(isParamValid("v_slip") ? &getFunctor<ADReal>("v_slip") : nullptr),
    _w_slip(isParamValid("w_slip") ? &getFunctor<ADReal>("w_slip") : nullptr),
    _add_slip_model(isParamValid("u_slip"))
{
  Moose::FV::setInterpolationMethod(
      *this, _advected_interp_method, "advected_interp_method");
}

Real
MSRLinearFVScalarAdvectionGenericRC::computeElemMatrixContribution()
{
  return _advected_interp_coeffs.first * _volumetric_face_flux * _current_face_area;
}

Real
MSRLinearFVScalarAdvectionGenericRC::computeNeighborMatrixContribution()
{
  return _advected_interp_coeffs.second * _volumetric_face_flux * _current_face_area;
}

Real
MSRLinearFVScalarAdvectionGenericRC::computeElemRightHandSideContribution()
{
  return 0.0;
}

Real
MSRLinearFVScalarAdvectionGenericRC::computeNeighborRightHandSideContribution()
{
  return 0.0;
}

Real
MSRLinearFVScalarAdvectionGenericRC::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "Expected a LinearFVAdvectionDiffusionBC.");

  const auto boundary_value_matrix_contrib =
      adv_bc->computeBoundaryValueMatrixContribution();

  const auto factor =
      (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  return boundary_value_matrix_contrib * factor * _volumetric_face_flux *
         _current_face_area;
}

Real
MSRLinearFVScalarAdvectionGenericRC::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "Expected a LinearFVAdvectionDiffusionBC.");

  const auto factor =
      (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  const auto boundary_value_rhs_contrib =
      adv_bc->computeBoundaryValueRHSContribution();

  return -boundary_value_rhs_contrib * factor * _volumetric_face_flux *
         _current_face_area;
}

void
MSRLinearFVScalarAdvectionGenericRC::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  _volumetric_face_flux =
      _mass_flux_provider.getVolumetricFaceFlux(Moose::FV::InterpMethod::RhieChow,
                                                *face_info,
                                                determineState(),
                                                _tid,
                                                false);

  if (_u_slip && face_info->neighborPtr())
  {
    const auto state = determineState();

    Moose::FaceArg face_arg{face_info,
                            Moose::FV::LimiterType::CentralDifference,
                            true,
                            false,
                            face_info->neighborPtr(),
                            nullptr};

    RealVectorValue velocity_slip_vel_vec;

    if (_u_slip)
      velocity_slip_vel_vec(0) = (*_u_slip)(face_arg, state).value();

    if (_v_slip)
      velocity_slip_vel_vec(1) = (*_v_slip)(face_arg, state).value();

    if (_w_slip)
      velocity_slip_vel_vec(2) = (*_w_slip)(face_arg, state).value();

    _volumetric_face_flux += velocity_slip_vel_vec * face_info->normal();
  }

  _advected_interp_coeffs =
      interpCoeffs(_advected_interp_method,
                   *_current_face_info,
                   true,
                   _volumetric_face_flux);
}
