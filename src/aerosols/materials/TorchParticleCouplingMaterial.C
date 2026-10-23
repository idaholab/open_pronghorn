//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#if defined(MOOSE_LIBTORCH_ENABLED) || defined(LIBTORCH_ENABLED)

#include "TorchParticleCouplingMaterial.h"

registerMooseObject("OpenPronghornApp", TorchParticleCouplingMaterial);

InputParameters
TorchParticleCouplingMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Material that exposes particle->fluid coupling (from TorchParticleCloudUserObject) "
      "as FV-accessible material properties and provides per-cell particle concentrations for visualization.");

  params.addRequiredParam<UserObjectName>(
      "particle_cloud_userobject",
      "Name of the TorchParticleCloudUserObject that computed the coupling fields.");

  params.addParam<bool>("output_force_density",
                        true,
                        "If true, output force density [N/m^3]. If false, output total cell force [N].");

  params.addParam<bool>("debug", false, "Print debug output for a few elements.");

  params.addParam<MooseFunctorName>("fx_name", "particle_fx", "Property name for x-force (or force density).");
  params.addParam<MooseFunctorName>("fy_name", "particle_fy", "Property name for y-force (or force density).");
  params.addParam<MooseFunctorName>("fz_name", "particle_fz", "Property name for z-force (or force density).");

  params.addParam<MooseFunctorName>("number_conc_name",
                                   "aerosol_number_conc",
                                   "Property name for aerosol number concentration [#/m^3].");
  params.addParam<MooseFunctorName>("mass_conc_name",
                                   "aerosol_mass_conc",
                                   "Property name for aerosol mass concentration [kg/m^3].");

  return params;
}

TorchParticleCouplingMaterial::TorchParticleCouplingMaterial(const InputParameters & parameters)
  : Material(parameters),
    _mesh_dimension(_mesh.dimension()),
    _output_force_density(getParam<bool>("output_force_density")),
    _debug(getParam<bool>("debug")),
    _particle_cloud(getUserObject<TorchParticleCloudUserObject>("particle_cloud_userobject")),
    _fx_name(getParam<MooseFunctorName>("fx_name")),
    _fy_name(getParam<MooseFunctorName>("fy_name")),
    _fz_name(getParam<MooseFunctorName>("fz_name")),
    _nconc_name(getParam<MooseFunctorName>("number_conc_name")),
    _mconc_name(getParam<MooseFunctorName>("mass_conc_name"))
{
  // Declare scalar properties (cell-wise constants in FV)
  _fx_prop = &declareGenericPropertyByName<Real, false>(_fx_name);
  _fy_prop = &declareGenericPropertyByName<Real, false>(_fy_name);
  _fz_prop = &declareGenericPropertyByName<Real, false>(_fz_name);

  _nconc_prop = &declareGenericPropertyByName<Real, false>(_nconc_name);
  _mconc_prop = &declareGenericPropertyByName<Real, false>(_mconc_name);
}

void
TorchParticleCouplingMaterial::initQpStatefulProperties()
{
  computeValues();
}

void
TorchParticleCouplingMaterial::computeQpProperties()
{
  computeValues();
}

void
TorchParticleCouplingMaterial::computeValues()
{
  if (!_current_elem)
  {
    (*_fx_prop)[_qp] = 0.0;
    (*_fy_prop)[_qp] = 0.0;
    (*_fz_prop)[_qp] = 0.0;
    (*_nconc_prop)[_qp] = 0.0;
    (*_mconc_prop)[_qp] = 0.0;
    return;
  }

  libMesh::VectorValue<Real> f(0.0, 0.0, 0.0);

  if (_output_force_density)
    f = _particle_cloud.cellForceDensity(_current_elem); // [N/m^3]
  else
    f = _particle_cloud.cellForce(_current_elem);        // [N]

  const Real fx = f(0);
  const Real fy = (_mesh_dimension > 1) ? f(1) : 0.0;
  const Real fz = (_mesh_dimension > 2) ? f(2) : 0.0;

  (*_fx_prop)[_qp] = fx;
  (*_fy_prop)[_qp] = fy;
  (*_fz_prop)[_qp] = fz;

  // Concentrations for visualization
  (*_nconc_prop)[_qp] = _particle_cloud.cellNumberConcentration(_current_elem);
  (*_mconc_prop)[_qp] = _particle_cloud.cellMassConcentration(_current_elem);

  if (_debug && _qp == 0)
  {
    _console << "[TorchParticleCouplingMaterial] elem=" << _current_elem->id()
             << "  f=(" << fx << ", " << fy << ", " << fz << ")"
             << (_output_force_density ? " [N/m^3]" : " [N]")
             << "  nconc=" << (*_nconc_prop)[_qp] << " [#/m^3]"
             << "  mconc=" << (*_mconc_prop)[_qp] << " [kg/m^3]"
             << std::endl;
  }
}

#endif
