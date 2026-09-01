//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVSRFMomentumBoussinesq.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "NS.h"
#include "FEProblemBase.h"
#include "SRFUtils.h"

registerMooseObject("OpenPronghornApp", LinearFVSRFMomentumBoussinesq);

InputParameters
LinearFVSRFMomentumBoussinesq::validParams()
{
  InputParameters params = LinearFVMomentumBoussinesq::validParams();
  params.addClassDescription("Represents the Boussinesq term in the SRF Navier Stokes momentum "
                             "equations, added to the right hand side.");
  // params.addParam<VariableName>(NS::T_fluid, "The fluid temperature variable.");
  // params.addRequiredParam<RealVectorValue>("gravity", "Gravitational acceleration vector.");
  // params.addParam<MooseFunctorName>("alpha_name",
  //                                   NS::alpha_boussinesq,
  //                                   "The name of the thermal expansion coefficient"
  //                                   "this is of the form rho = rho_ref*(1-alpha (T-T_ref))");
  // params.addRequiredParam<Real>("ref_temperature", "The value for the reference temperature.");
  // params.addRequiredParam<MooseFunctorName>(NS::density, "The value for the density");
  // MooseEnum momentum_component("x=0 y=1 z=2");
  // params.addRequiredParam<MooseEnum>(
  //     "momentum_component",
  //     momentum_component,
  //     "The component of the momentum equation that this kernel applies to.");
  params.addRequiredParam<MooseFunctorName>("pitch_angle", "The pitch angle of rotation between the metacenter and the body reference frame.");
  params.addRequiredParam<MooseFunctorName>("yaw_angle", "The yaw angle of rotation between the metacenter and the body reference frame.");
  params.addRequiredParam<MooseFunctorName>("roll_angle", "The roll angle of rotation between the metacenter and the body reference frame.");
  return params;
}

LinearFVSRFMomentumBoussinesq::LinearFVSRFMomentumBoussinesq(const InputParameters & params)
  : LinearFVMomentumBoussinesq(params),
    _pitch_angle(getFunctor<Real>("pitch_angle")),
    _yaw_angle(getFunctor<Real>("yaw_angle")),
    _roll_angle(getFunctor<Real>("roll_angle"))
{
}

Real
LinearFVSRFMomentumBoussinesq::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVSRFMomentumBoussinesq::computeRightHandSideContribution()
{
  const auto elem = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();

  const RealVectorValue gravity_brf = NS::SRF::rotateVectorInertialToBody(_gravity,
                                           _pitch_angle(elem,state), _yaw_angle(elem,state), _roll_angle(elem,state));

  return -_alpha(elem, state) * gravity_brf(_index) * _rho(elem, state) *
         (_temperature_var.getElemValue(*_current_elem_info, state) - _ref_temperature) *
         _current_elem_volume;
}
