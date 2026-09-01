//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVSRFMomentumHeave.h"
#include "NS.h"

#include <cmath>

registerMooseObject("OpenPronghornApp", LinearFVSRFMomentumHeave);

InputParameters
LinearFVSRFMomentumHeave::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();

  params.addClassDescription(
      "Adds rho times the fictitious acceleration associated with sinusoidal heave to a selected "
      "body-frame momentum equation.");

  params.addRequiredParam<MooseFunctorName>(
      "heave_amp", "Heave displacement amplitude in the selected body-frame direction [m].");
  params.addRequiredParam<MooseFunctorName>("heave_per", "Heave period [s].");
  params.addParam<MooseFunctorName>(
      "heave_pha", 0.0, "Heave displacement phase shift [degrees].");

  params.addRequiredParam<MooseFunctorName>(NS::density, "Fluid density [kg/m^3].");

  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "Body-frame momentum component to which the heave source is applied.");

  return params;
}

LinearFVSRFMomentumHeave::LinearFVSRFMomentumHeave(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _index(getParam<MooseEnum>("momentum_component")),
    _heave_amp(getFunctor<Real>("heave_amp")),
    _heave_per(getFunctor<Real>("heave_per")),
    _heave_pha(getFunctor<Real>("heave_pha")),
    _rho(getFunctor<Real>(NS::density))
{
}

Real
LinearFVSRFMomentumHeave::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVSRFMomentumHeave::computeRightHandSideContribution()
{
  const auto elem = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();

  const Real amplitude = _heave_amp(elem, state);

  if (!std::isfinite(amplitude))
    mooseError(name(), ": heave_amp must be finite.");

  // A zero amplitude disables heave and permits a zero placeholder period.
  if (amplitude == 0.0)
    return 0.0;

  const Real period = _heave_per(elem, state);
  if (!std::isfinite(period) || period <= 0.0)
    mooseError(name(),
               ": heave_per must be finite and greater than zero when heave_amp is nonzero. "
               "Received ",
               period,
               ".");

  const Real phase_deg = _heave_pha(elem, state);
  if (!std::isfinite(phase_deg))
    mooseError(name(), ": heave_pha must be finite.");

  const Real omega = 2.0 * libMesh::pi / period;
  const Real phase_rad = phase_deg * libMesh::pi / 180.0;

  // eta = A sin(omega t + phi)
  // a_frame = eta_ddot = -A omega^2 sin(omega t + phi)
  // a_fluid = -a_frame = +A omega^2 sin(omega t + phi)
  const Real fictitious_heave_acceleration =
      amplitude * omega * omega * std::sin(omega * _t + phase_rad);

  // The scalar heave motion is applied along the selected body-frame component.
  RealVectorValue heave_acceleration_brf(0.0, 0.0, 0.0);
  heave_acceleration_brf(_index) = fictitious_heave_acceleration;

  return _rho(elem, state) * heave_acceleration_brf(_index) * _current_elem_volume;
}
