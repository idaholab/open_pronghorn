#include "MSRLinearRhieChowMassFlux.h"
#include "LinearSystem.h"

#include "libmesh/dof_map.h"
#include "libmesh/linear_implicit_system.h"

#include <cmath>
#include <limits>
#include <vector>

registerMooseObject("OpenPronghornApp", MSRLinearRhieChowMassFlux);

InputParameters
MSRLinearRhieChowMassFlux::validParams()
{
  auto params = RhieChowMassFlux::validParams();
  params.addClassDescription(
      "Linear SIMPLE Rhie-Chow mass-flux provider with read-only access to the effective "
      "cell momentum diagonal for MultiApp Rhie-Chow reconstruction.");
  return params;
}

MSRLinearRhieChowMassFlux::MSRLinearRhieChowMassFlux(const InputParameters & params)
  : RhieChowMassFlux(params)
{
}

Real
MSRLinearRhieChowMassFlux::momentumDiagonal(const Elem & elem, const unsigned int component) const
{
  if (component >= _dim)
    mooseError(name(), ": requested momentum component ", component, " but mesh dimension is ", _dim);

  if (component >= _Ainv_raw.size() || !_Ainv_raw[component])
    mooseError(name(),
               ": Rhie-Chow cell coefficient data are not available yet. "
               "Evaluate the bridge after the SIMPLE momentum/pressure solve (for example at "
               "TIMESTEP_END).");

  if (!_cell_volumes || !_pressure_system)
    mooseError(name(), ": Rhie-Chow mesh/system information has not been initialized.");

  std::vector<dof_id_type> momentum_dofs;
  _momentum_implicit_systems[component]->get_dof_map().dof_indices(
      &elem, momentum_dofs, _vel[component]->number());

  if (momentum_dofs.size() != 1)
    mooseError(name(),
               ": expected exactly one FV momentum DOF for element ",
               elem.id(),
               " and component ",
               component,
               ", found ",
               momentum_dofs.size(),
               ".");

  std::vector<dof_id_type> pressure_dofs;
  _pressure_system->system().get_dof_map().dof_indices(&elem, pressure_dofs, _p->number());

  if (pressure_dofs.size() != 1)
    mooseError(name(),
               ": expected exactly one FV pressure DOF for element ",
               elem.id(),
               ", found ",
               pressure_dofs.size(),
               ".");

  const Real volume_over_A = (*_Ainv_raw[component])(momentum_dofs[0]);
  const Real cell_volume = (*_cell_volumes)(pressure_dofs[0]);

  if (!std::isfinite(volume_over_A) || !std::isfinite(cell_volume))
    mooseError(name(), ": non-finite Rhie-Chow coefficient data on element ", elem.id(), ".");

  if (std::abs(volume_over_A) <= std::numeric_limits<Real>::epsilon())
    mooseError(name(),
               ": cannot recover the momentum diagonal on element ",
               elem.id(),
               " because V/A is zero (or numerically indistinguishable from zero).");

  // RhieChowMassFlux stores Ainv_raw = V/A after computeHbyA().
  // INSFVRhieChowInterpolator expects a = A and reconstructs D = V/a.
  return cell_volume / volume_over_A;
}
