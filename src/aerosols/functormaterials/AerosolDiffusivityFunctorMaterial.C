//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AerosolDiffusivityFunctorMaterial.h"
#include "MooseMesh.h"
#include "NS.h"

#include "libmesh/libmesh_common.h" // for libMesh::pi

registerMooseObject("OpenPronghornApp", AerosolDiffusivityFunctorMaterial);

InputParameters
AerosolDiffusivityFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Computes Brownian and turbulent aerosol diffusion coefficients as functors.");

  // Fluid properties
  params.addRequiredParam<MooseFunctorName>(NS::density, "Carrier fluid density functor [kg/m^3].");
  params.addRequiredParam<MooseFunctorName>(NS::mu,
                                            "Carrier fluid dynamic viscosity functor [Pa·s].");
  params.addRequiredParam<MooseFunctorName>(NS::T_fluid, "Carrier fluid temperature functor [K].");

  // Turbulent viscosity (assumed provided by another material, e.g. k-epsilon)
  params.addRequiredParam<MooseFunctorName>(NS::mu_t,
                                            "Turbulent dynamic viscosity functor [Pa·s].");

  // Particle properties
  params.addRequiredParam<Real>("particle_diameter", "Aerosol particle diameter [m].");
  params.addParam<Real>("mean_free_path",
                        0.066e-6,
                        "Mean free path of the carrier gas [m] used in the "
                        "Cunningham slip correction.");

  // Turbulent Schmidt number
  params.addParam<Real>("Sc_t", 1.0, "Turbulent Schmidt number Sc_t for the aerosol.");

  return params;
}

AerosolDiffusivityFunctorMaterial::AerosolDiffusivityFunctorMaterial(
    const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _T(getFunctor<Real>(NS::T_fluid)),
    _mu_t(getFunctor<Real>(NS::mu_t)),
    _d_p(getParam<Real>("particle_diameter")),
    _lambda_mfp(getParam<Real>("mean_free_path")),
    _Sc_t(getParam<Real>("Sc_t"))
{
  // Low limits to avoid zero diffusion (can be tuned as needed)
  const Real D_B_low_limit = 1e-20;
  const Real D_t_low_limit = 1e-20;

  // Boltzmann constant [J/K]
  const Real k_B = 1.380649e-23;

  // Brownian diffusivity functor
  addFunctorProperty<Real>(
      "D_B",
      [this, D_B_low_limit, k_B](const auto & r, const auto & t) -> Real
      {
        using std::exp;
        using std::max;

        const Real mu = _mu(r, t);
        const Real T = _T(r, t);

        if (_d_p <= 0.0 || _lambda_mfp <= 0.0 || mu <= 0.0 || T <= 0.0)
          return D_B_low_limit;

        // Cunningham slip correction factor
        const Real Kn_ratio = _lambda_mfp / _d_p;
        const Real C_c =
            Real(1.0) +
            Kn_ratio * (Real(2.34) + Real(1.05) * exp(Real(-0.39) * _d_p / _lambda_mfp));

        const Real D_B_raw = k_B * T * C_c / (Real(3.0) * libMesh::pi * mu * _d_p);

        return max(D_B_low_limit, D_B_raw);
      });

  // Turbulent diffusivity functor
  addFunctorProperty<Real>("D_t",
                           [this, D_t_low_limit](const auto & r, const auto & t) -> Real
                           {
                             using std::max;

                             const Real rho = _rho(r, t);
                             const Real mu_t = _mu_t(r, t);

                             if (rho <= 0.0 || _Sc_t <= 0.0)
                               return D_t_low_limit;

                             const Real D_t_raw = mu_t / (rho * _Sc_t);

                             return max(D_t_low_limit, D_t_raw);
                           });

  // Effective diffusivity (Brownian + turbulent)
  addFunctorProperty<Real>("D_eff",
                           [this](const auto & r, const auto & t) -> Real
                           {
                             const Real DB = this->getFunctor<Real>("D_B")(r, t);
                             const Real Dt = this->getFunctor<Real>("D_t")(r, t);
                             return DB + Dt;
                           });
}
