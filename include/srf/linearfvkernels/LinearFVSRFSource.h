//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"
#include "NS.h"

/**
 * Kernel that adds contributions from a external source term discretized using the finite volume
 * method to a linear system.
 */
class LinearFVSRFSource : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVSRFSource(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:
  // Source-term vector
  const Moose::Functor<RealVectorValue> & _source_density_vector;
    // Scaling Factor.
  const Moose::Functor<Real> & _scale;
  // Direction (x,y,z) of the source component in the inertial reference frame.
  const MooseEnum _index;
  // Pitch angle
  const Moose::Functor<Real> & _pitch_angle;
  // Yaw angle
  const Moose::Functor<Real> & _yaw_angle;
  // Roll angle
  const Moose::Functor<Real> & _roll_angle;

};
