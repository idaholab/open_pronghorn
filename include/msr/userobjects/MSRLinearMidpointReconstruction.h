#pragma once

#include "GeneralUserObject.h"

/**
 * Reconstruct physical end-of-step states for all solver systems using
 * MSRLinearImplicitMidpoint.
 *
 * The custom linear midpoint solve stores the midpoint stage Y. For the
 * operator-split EPI3V chain the endpoint y_{n+1}=2Y-y_n must exist before
 * TIMESTEP_END transfers launch the chemistry MultiApp. MOOSE calls
 * TimeIntegrator::postStep() too late for that nested-MultiApp ordering, so
 * this object performs the reconstruction during EXEC_TIMESTEP_END.
 */
class MSRLinearMidpointReconstruction : public GeneralUserObject
{
public:
  static InputParameters validParams();
  MSRLinearMidpointReconstruction(const InputParameters & parameters);

  void initialize() override {}
  void execute() override;
  void finalize() override {}
};
