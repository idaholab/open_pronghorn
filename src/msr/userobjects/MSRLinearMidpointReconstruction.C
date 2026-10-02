#include "MSRLinearMidpointReconstruction.h"

#include "MSRLinearImplicitMidpoint.h"
#include "FEProblemBase.h"
#include "SolverSystem.h"

registerMooseObject("OpenPronghornApp", MSRLinearMidpointReconstruction);

InputParameters
MSRLinearMidpointReconstruction::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Reconstructs MSRLinearImplicitMidpoint endpoint states at TIMESTEP_END "
      "before MultiApp transfers.");
  return params;
}

MSRLinearMidpointReconstruction::MSRLinearMidpointReconstruction(
    const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}

void
MSRLinearMidpointReconstruction::execute()
{
  unsigned int reconstructed = 0;

  for (const auto i : make_range(_fe_problem.numSolverSystems()))
  {
    auto & sys = _fe_problem.getSolverSystem(i);

    for (auto & ti : sys.getTimeIntegrators())
      if (auto * midpoint = dynamic_cast<MSRLinearImplicitMidpoint *>(ti.get()))
      {
        midpoint->reconstructSolution();
        ++reconstructed;
      }
  }

  if (reconstructed == 0)
    mooseError(
        "MSRLinearMidpointReconstruction executed but found no "
        "MSRLinearImplicitMidpoint solver systems.");
}
