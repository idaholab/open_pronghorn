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

/**
 * Adds the inertial source associated with sinusoidal heave to one body-frame momentum equation.
 *
 * The prescribed body-frame displacement is
 *
 *   eta(t) = A sin(2 pi t / T + phi),
 *
 * where A is the displacement amplitude, T is the period, and phi is the phase. The acceleration
 * of the non-inertial frame is eta_ddot. The fictitious acceleration acting on the fluid is
 * -eta_ddot, so this kernel adds
 *
 *   rho A (2 pi / T)^2 sin(2 pi t / T + phi)
 *
 * to the selected momentum-equation right-hand side.
 */
class LinearFVSRFMomentumHeave : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  LinearFVSRFMomentumHeave(const InputParameters & params);

  virtual Real computeMatrixContribution() override;
  virtual Real computeRightHandSideContribution() override;

protected:
  /// Index of the body-frame momentum component: x=0, y=1, or z=2.
  const unsigned int _index;

  /// Heave displacement amplitude [m].
  const Moose::Functor<Real> & _heave_amp;

  /// Heave period [s].
  const Moose::Functor<Real> & _heave_per;

  /// Heave displacement phase [degrees].
  const Moose::Functor<Real> & _heave_pha;

  /// Fluid density [kg/m^3].
  const Moose::Functor<Real> & _rho;
};
