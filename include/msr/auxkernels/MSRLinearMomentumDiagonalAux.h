#pragma once

#include "AuxKernel.h"

class MSRLinearRhieChowMassFlux;

/**
 * Writes one component of the linear SIMPLE momentum coefficient A into an
 * elemental/FV auxiliary variable so it can be transferred to a MultiApp.
 */
class MSRLinearMomentumDiagonalAux : public AuxKernel
{
public:
  static InputParameters validParams();
  MSRLinearMomentumDiagonalAux(const InputParameters & params);

protected:
  virtual Real computeValue() override;

  const MSRLinearRhieChowMassFlux & _rc;
  const unsigned int _component;
};
