//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

/**
 * Converts temperature to enthalpy or enthalpy to temperature using
 * functor material properties. The derivatives are discarded, so can't
 * be used with AD.
 */
class LinearFVSRFFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();
  LinearFVSRFFunctorMaterial(const InputParameters & parameters);

protected:
  /// Metacenter origin
  const RealVectorValue _mc_origin;

  /// SRF input mode: fixed or pitch_yaw_roll
  const MooseEnum _SRF_input_mode;

  /// Fixed angles
  const Moose::Functor<Real> * _pitch_angle_fixed;
  const Moose::Functor<Real> * _yaw_angle_fixed;
  const Moose::Functor<Real> * _roll_angle_fixed;

  /// Fixed angular velocities
  const Moose::Functor<Real> & _pitch_omega_fixed;
  const Moose::Functor<Real> & _yaw_omega_fixed;
  const Moose::Functor<Real> & _roll_omega_fixed;

  /// Fixed angular accelerations
  const Moose::Functor<Real> & _pitch_omegadot_fixed;
  const Moose::Functor<Real> & _yaw_omegadot_fixed;
  const Moose::Functor<Real> & _roll_omegadot_fixed;

  /// Dynamic pitch-yaw-roll inputs
  const Moose::Functor<Real> * _pitch_amp;
  const Moose::Functor<Real> * _pitch_per;
  const Moose::Functor<Real> & _pitch_pha;

  const Moose::Functor<Real> * _yaw_amp;
  const Moose::Functor<Real> * _yaw_per;
  const Moose::Functor<Real> & _yaw_pha;

  const Moose::Functor<Real> * _roll_amp;
  const Moose::Functor<Real> * _roll_per;
  const Moose::Functor<Real> & _roll_pha;

  using UserObjectInterface::getUserObject;
};
