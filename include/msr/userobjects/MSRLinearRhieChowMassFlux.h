#pragma once

#include "RhieChowMassFlux.h"

/**
 * Small OpenPronghorn bridge around MOOSE's linear SIMPLE Rhie-Chow object.
 *
 * RhieChowMassFlux stores the cell pressure-correction coefficient as
 *     Ainv_raw = V / A
 * after computeHbyA(), where A is the effective momentum diagonal (or the
 * SIMPLEC row-sum equivalent) and V includes the coordinate factor.
 *
 * INSFVRhieChowInterpolator expects A itself when a_u/a_v are supplied from
 * auxiliary fields.  This subclass exposes A without changing the linear
 * SIMPLE/Rhie-Chow algorithm or the stored face mass flux.
 */
class MSRLinearRhieChowMassFlux : public RhieChowMassFlux
{
public:
  static InputParameters validParams();
  MSRLinearRhieChowMassFlux(const InputParameters & params);

  /// Return the effective cell momentum coefficient A for one velocity component.
  Real momentumDiagonal(const Elem & elem, unsigned int component) const;
};
