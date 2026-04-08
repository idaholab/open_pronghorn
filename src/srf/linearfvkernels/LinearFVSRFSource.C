//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVSRFSource.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "SRFUtils.h"

registerMooseObject("OpenPronghornApp", LinearFVSRFSource);

InputParameters
LinearFVSRFSource::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Represents the transformed (solution-independent) momentum source term vector "
      "from the metacenter (inertial) reference frame to the body reference frame.");

  params.addRequiredParam<MooseFunctorName>("source_density_vector", "The source density vector in the metacenter (inertial) reference frame.");
  params.addRequiredParam<MooseFunctorName>("scaling_factor", "Coefficient to multiply the body force term with.");
    MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>("momentum_component",momentum_component,"The component in the body reference frame of the momentum equation that this kernel applies to.");
  params.addRequiredParam<MooseFunctorName>("pitch_angle", "The pitch angle of rotation between the metacenter and the body reference frame.");
  params.addRequiredParam<MooseFunctorName>("yaw_angle", "The yaw angle of rotation between the metacenter and the body reference frame.");
  params.addRequiredParam<MooseFunctorName>("roll_angle", "The roll angle of rotation between the metacenter and the body reference frame.");

  return params;
}

LinearFVSRFSource::LinearFVSRFSource(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _source_density_vector(getFunctor<RealVectorValue>("source_density_vector")),
    _scale(getFunctor<Real>("scaling_factor")),
    _index(getParam<MooseEnum>("momentum_component")),
    _pitch_angle(getFunctor<Real>("pitch_angle")),
    _yaw_angle(getFunctor<Real>("yaw_angle")),
    _roll_angle(getFunctor<Real>("roll_angle"))
{
}

Real
LinearFVSRFSource::computeMatrixContribution()
{
  // This doesn't contribute to the matrix as we assumed it was independent of
  // the solution
  return 0.0;
}

Real
LinearFVSRFSource::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto state_arg = determineState();

  const RealVectorValue source_density_vector_brf = NS::SRF::rotateVectorInertialToBody(_source_density_vector(elem_arg,state_arg),
                                           _pitch_angle(elem_arg,state_arg), _yaw_angle(elem_arg,state_arg), _roll_angle(elem_arg,state_arg));

  return _scale(elem_arg, state_arg) * source_density_vector_brf(_index) * _current_elem_volume;
}
