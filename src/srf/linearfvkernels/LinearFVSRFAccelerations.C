//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVSRFAccelerations.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "SRFUtils.h"
#include "FEProblemBase.h"

registerMooseObject("OpenPronghornApp", LinearFVSRFAccelerations);

InputParameters
LinearFVSRFAccelerations::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Represents the accelerations on the body reference frame "
      "from the metacenter reference frame motion.");

  params.addRequiredParam<MooseFunctorName>(NS::density, "Density");
  params.addRequiredParam<MooseFunctorName>("omega_brf", "The angular velocity in the body reference frame.");
  params.addRequiredParam<MooseFunctorName>("omega_dot_brf", "The angular acceleration in the body reference frame.");
  params.addRequiredParam<MooseFunctorName>("r_mc", "Vector Coordinates of the metacenter reference frame w.r.t the body origin.");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>("momentum_component",momentum_component,"The component in the body reference frame of the momentum equation that this kernel applies to.");

  params.addRequiredParam<SolverVariableName>("u", "The velocity in the x direction.");
  params.addRequiredParam<SolverVariableName>("v", "The velocity in the y direction.");
  params.addParam<SolverVariableName>("w", "The velocity in the z direction.");
  return params;
}

LinearFVSRFAccelerations::LinearFVSRFAccelerations(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _rho(getFunctor<Real>(NS::density)),
    _omega_brf(getFunctor<RealVectorValue>("omega_brf")),
    _omega_dot_brf(getFunctor<RealVectorValue>("omega_dot_brf")),
    _index(getParam<MooseEnum>("momentum_component")),
    _r_mc(getFunctor<RealVectorValue>("r_mc")),
    _u_var(dynamic_cast<const MooseLinearVariableFVReal *>(
        &_fe_problem.getVariable(_tid, getParam<SolverVariableName>("u")))),
    _v_var(dynamic_cast<const MooseLinearVariableFVReal *>(
        &_fe_problem.getVariable(_tid, getParam<SolverVariableName>("v")))),
    _w_var(params.isParamValid("w")
               ? dynamic_cast<const MooseLinearVariableFVReal *>(
                     &_fe_problem.getVariable(_tid, getParam<SolverVariableName>("w"))): nullptr)
{
}

Real
LinearFVSRFAccelerations::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVSRFAccelerations::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto state_arg = determineState();

  const Real rho = _rho(elem_arg, state_arg);

  // Centripetal acceleration
  const RealVectorValue omega = _omega_brf(elem_arg, state_arg);
  const RealVectorValue rmc = _r_mc(elem_arg, state_arg);
  const RealVectorValue omega_x_omega_x_r = omega.cross(omega.cross(rmc));
  // Tangential acceleration
  const RealVectorValue omega_dot = _omega_dot_brf(elem_arg, state_arg);
  const RealVectorValue omega_dot_x_r = omega_dot.cross(rmc);
  // Coriolis effect
  const RealVectorValue vel((*_u_var).getElemValue(*_current_elem_info, state_arg),
                            (*_v_var).getElemValue(*_current_elem_info, state_arg),
                            _w_var ? (*_w_var).getElemValue(*_current_elem_info, state_arg) : 0.0);

  const RealVectorValue coriolis = 2.0 * omega.cross(vel);

  return -rho * (omega_x_omega_x_r(_index) + omega_dot_x_r(_index) + coriolis(_index)) * _current_elem_volume;
}
