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
class LinearFVSRFAccelerations : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVSRFAccelerations(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:
  /// Density
  const Moose::Functor<Real> & _rho;
  // Body reference frame angular velocity vector
  const Moose::Functor<RealVectorValue> & _omega_brf;
  // Body reference frame angular acceleration vector
  const Moose::Functor<RealVectorValue> & _omega_dot_brf;
  // Direction (x,y,z) of the source component in the inertial reference frame.
  const MooseEnum _index;
  // Metacenter origin vector
  const Moose::Functor<RealVectorValue> & _r_mc;
  /// x-velocity
  const Moose::Functor<Real> & _u_var;
  /// y-velocity
  const Moose::Functor<Real> & _v_var;
  /// z-velocity
  const Moose::Functor<Real> & _w_var;
  // Flag to add Coriolis or not.
  const bool _coriolis;

};
