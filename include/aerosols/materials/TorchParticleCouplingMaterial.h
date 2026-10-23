//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#if defined(MOOSE_LIBTORCH_ENABLED) || defined(LIBTORCH_ENABLED)

#include "Material.h"
#include "MooseTypes.h"
#include "TorchParticleCloudUserObject.h"

/**
 * TorchParticleCouplingMaterial
 *
 * Exposes particle->fluid coupling and particle-derived visualization fields
 * from TorchParticleCloudUserObject as Material properties (cell-wise constants).
 *
 * Properties (default names):
 *  - particle_fx, particle_fy, particle_fz : force density [N/m^3] (or force [N])
 *  - aerosol_number_conc : [#/m^3]
 *  - aerosol_mass_conc   : [kg/m^3]
 *
 * This material does not do any MPI handling; it reads whatever the local userobject provides.
 */
class TorchParticleCouplingMaterial : public Material
{
public:
  static InputParameters validParams();
  TorchParticleCouplingMaterial(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

private:
  void computeValues();

  const unsigned int _mesh_dimension;

  const bool _output_force_density;
  const bool _debug;

  const TorchParticleCloudUserObject & _particle_cloud;

  const MooseFunctorName & _fx_name;
  const MooseFunctorName & _fy_name;
  const MooseFunctorName & _fz_name;

  const MooseFunctorName & _nconc_name;
  const MooseFunctorName & _mconc_name;

  GenericMaterialProperty<Real, false> * _fx_prop = nullptr;
  GenericMaterialProperty<Real, false> * _fy_prop = nullptr;
  GenericMaterialProperty<Real, false> * _fz_prop = nullptr;

  GenericMaterialProperty<Real, false> * _nconc_prop = nullptr;
  GenericMaterialProperty<Real, false> * _mconc_prop = nullptr;
};

#endif
