//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVMomentumBoussinesq.h"


/**
 * Kernel that adds contributions from a external source term discretized using the finite volume
 * method to a linear system.
 */
class LinearFVSRFMomentumBoussinesq : public LinearFVMomentumBoussinesq
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVSRFMomentumBoussinesq(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:
  // Pitch angle
  const Moose::Functor<Real> & _pitch_angle;
  // Yaw angle
  const Moose::Functor<Real> & _yaw_angle;
  // Roll angle
  const Moose::Functor<Real> & _roll_angle;

};
