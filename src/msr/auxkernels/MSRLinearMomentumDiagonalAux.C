#include "MSRLinearMomentumDiagonalAux.h"
#include "MSRLinearRhieChowMassFlux.h"

registerMooseObject("OpenPronghornApp", MSRLinearMomentumDiagonalAux);

InputParameters
MSRLinearMomentumDiagonalAux::validParams()
{
  auto params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object", "MSRLinearRhieChowMassFlux object from the linear SIMPLE solve.");
  params.addRequiredParam<MooseEnum>(
      "component", MooseEnum("x y z"), "Momentum component whose effective diagonal is exported.");
  params.addClassDescription(
      "Exports a linear SIMPLE momentum diagonal coefficient for MultiApp Rhie-Chow reconstruction.");
  return params;
}

MSRLinearMomentumDiagonalAux::MSRLinearMomentumDiagonalAux(const InputParameters & params)
  : AuxKernel(params),
    _rc(getUserObject<MSRLinearRhieChowMassFlux>("rhie_chow_user_object")),
    _component(getParam<MooseEnum>("component") == "x"
                   ? 0
                   : (getParam<MooseEnum>("component") == "y" ? 1 : 2))
{
}

Real
MSRLinearMomentumDiagonalAux::computeValue()
{
  return _rc.momentumDiagonal(*_current_elem, _component);
}
