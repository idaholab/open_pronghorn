//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#include "AdvancedKEpsilonPhysics.h"

#include "WCNSFVFlowPhysicsBase.h"
#include "WCNSFVFluidHeatTransferPhysicsBase.h"
#include "kEpsilonViscosity.h"
#include "kEpsilonTKESourceSink.h"
#include "kEpsilonTKEDSourceSink.h"
#include "WallDistanceAux.h"

registerWCNSFVTurbulenceBaseTasks("OpenPronghornApp", AdvancedKEpsilonPhysics);
registerMooseAction("OpenPronghornApp",
                    AdvancedKEpsilonPhysics,
                    "add_interpolation_method_physics");
registerMooseAction("OpenPronghornApp", AdvancedKEpsilonPhysics, "add_functor_material");

// Parameters with the same name and meaning in the three OpenPronghorn k-epsilon objects, and
// which can therefore be forwarded verbatim. Each object only picks up the subset it declares.
const std::vector<std::string> AdvancedKEpsilonPhysics::_variant_parameters = {
    "k_epsilon_variant",
    "two_layer_flavor",
    "scale_limiter",
    "nonlinear_model",
    "curvature_model",
    "use_buoyancy",
    "use_compressibility",
    "use_curvature_correction",
    "use_yap",
    "use_low_re_Gprime",
    "Cd0",
    "Cd1",
    "Cd2",
    "Ca0",
    "Ca1",
    "Ca2",
    "Ca3",
    "C3_eps",
    "C_lowRe",
    "D_lowRe",
    "E_lowRe",
    "Cw",
    "C_M",
    "speed_of_sound"};

InputParameters
AdvancedKEpsilonPhysics::validParams()
{
  InputParameters params = WCNSLinearFVTurbulencePhysics::validParams();
  params.addClassDescription(
      "Define a k-epsilon turbulence model for an incompressible or weakly-compressible Navier "
      "Stokes flow with a linear finite volume discretization, including the two-layer, low "
      "Reynolds number and realizable variants and their associated corrections");

  const auto viscosity_params = kEpsilonViscosity::validParams();
  const auto tke_params = kEpsilonTKESourceSink::validParams();
  const auto tked_params = kEpsilonTKEDSourceSink::validParams();

  // Model variant selection and near-wall formulation
  params.transferParam<MooseEnum>(viscosity_params, "k_epsilon_variant");
  params.transferParam<MooseEnum>(viscosity_params, "two_layer_flavor");
  params.transferParam<MooseEnum>(viscosity_params, "scale_limiter");
  params.transferParam<MooseEnum>(tke_params, "nonlinear_model");
  params.transferParam<MooseEnum>(tke_params, "curvature_model");

  // Optional corrections to the k and epsilon equations
  params.transferParam<bool>(tke_params, "use_buoyancy");
  params.transferParam<bool>(tke_params, "use_compressibility");
  params.transferParam<bool>(tke_params, "use_curvature_correction");
  params.transferParam<bool>(tked_params, "use_yap");
  params.transferParam<bool>(tked_params, "use_low_re_Gprime");

  // Closure coefficients of the variants
  params.transferParam<Real>(viscosity_params, "Cd0");
  params.transferParam<Real>(viscosity_params, "Cd1");
  params.transferParam<Real>(viscosity_params, "Cd2");
  params.transferParam<Real>(viscosity_params, "Ca0");
  params.transferParam<Real>(viscosity_params, "Ca1");
  params.transferParam<Real>(viscosity_params, "Ca2");
  params.transferParam<Real>(viscosity_params, "Ca3");
  params.transferParam<Real>(tked_params, "C3_eps");
  params.transferParam<Real>(tked_params, "C_lowRe");
  params.transferParam<Real>(tked_params, "D_lowRe");
  params.transferParam<Real>(tked_params, "E_lowRe");
  params.transferParam<Real>(tked_params, "Cw");
  params.transferParam<Real>(tke_params, "C_M");
  params.transferParam<MooseFunctorName>(tke_params, "speed_of_sound");

  // 'Ct' has a different meaning in the turbulent viscosity and in the epsilon equation, so the
  // two cannot share a single parameter
  params.transferParam<Real>(
      viscosity_params,
      "Ct",
      "time_scale_Ct",
      "Time-scale constant Ct used in the k-epsilon time scale T = max(Te, Ct*sqrt(nu/epsilon)) "
      "when scale_limiter='standard'.");
  params.transferParam<Real>(tked_params,
                             "Ct",
                             "yap_Ct",
                             "Ct coefficient used in the Yap / ambient epsilon source terms.");

  // The two-layer, low-Reynolds and Yap formulations need the distance to the nearest wall. It
  // is created by this Physics unless the user points to an existing functor
  params.addParam<MooseFunctorName>(
      NS::wall_distance,
      NS::wall_distance,
      "Functor providing the distance to the nearest wall, used by the two-layer, low Reynolds "
      "number and Yap formulations. If not specified, an auxiliary variable of this name is "
      "created and computed on the turbulence walls");

  params.addParamNamesToGroup("k_epsilon_variant two_layer_flavor nonlinear_model curvature_model "
                              "use_buoyancy use_compressibility use_curvature_correction use_yap "
                              "use_low_re_Gprime",
                              "K-Epsilon model");
  params.addParamNamesToGroup("Cd0 Cd1 Cd2 Ca0 Ca1 Ca2 Ca3 C3_eps C_lowRe D_lowRe E_lowRe Cw C_M "
                              "time_scale_Ct yap_Ct speed_of_sound scale_limiter wall_distance",
                              "K-Epsilon model");

  return params;
}

AdvancedKEpsilonPhysics::AdvancedKEpsilonPhysics(const InputParameters & parameters)
  : WCNSLinearFVTurbulencePhysics(parameters),
    _wall_distance_name(getParam<MooseFunctorName>(NS::wall_distance)),
    _create_wall_distance(needsWallDistance() && !isParamSetByUser(NS::wall_distance))
{
  if (_turbulence_model != "k-epsilon")
    errorDependentParameter("turbulence_handling", "k-epsilon", _variant_parameters);

  if (_create_wall_distance)
    saveAuxVariableName(_wall_distance_name);
}

bool
AdvancedKEpsilonPhysics::needsWallDistance() const
{
  const auto variant = getParam<MooseEnum>("k_epsilon_variant");
  return variant == "StandardTwoLayer" || variant == "RealizableTwoLayer" ||
         variant == "StandardLowRe" || getParam<bool>("use_yap") ||
         getParam<bool>("use_low_re_Gprime");
}

Real
AdvancedKEpsilonPhysics::getConstantCoefficient(const std::string & param_name) const
{
  const auto & functor_name = getParam<MooseFunctorName>(param_name);
  if (!MooseUtils::isFloat(functor_name))
    paramError(param_name,
               "The OpenPronghorn k-epsilon objects only accept a constant value for this "
               "coefficient, but '",
               functor_name,
               "' is not a constant");
  return MooseUtils::convert<Real>(functor_name);
}

void
AdvancedKEpsilonPhysics::addAuxiliaryVariables()
{
  WCNSLinearFVTurbulencePhysics::addAuxiliaryVariables();

  if (!_create_wall_distance)
    return;

  // WallDistanceAux uses the nonlinear finite volume machinery, so the wall distance cannot be a
  // MooseLinearVariableFVReal like the other auxiliary variables of this Physics
  const std::string var_type = "MooseVariableFVReal";
  auto params = getFactory().getValidParams(var_type);
  assignBlocks(params, _blocks);
  if (shouldCreateVariable(_wall_distance_name, _blocks, /*error if aux*/ false))
    getProblem().addAuxVariable(var_type, _wall_distance_name, params);
}

void
AdvancedKEpsilonPhysics::addAuxiliaryKernels()
{
  WCNSLinearFVTurbulencePhysics::addAuxiliaryKernels();

  if (!_create_wall_distance)
    return;

  const std::string kernel_type = "WallDistanceAux";
  auto params = getFactory().getValidParams(kernel_type);
  assignBlocks(params, _blocks);
  params.set<AuxVariableName>("variable") = _wall_distance_name;
  params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
  // The distance to the walls only depends on the mesh
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL};
  getProblem().addAuxKernel(kernel_type, prefix() + "wall_distance_aux", params);
}

void
AdvancedKEpsilonPhysics::addTurbulentViscosityAuxKernel()
{
  const std::string u_names[3] = {"u", "v", "w"};
  const std::string kernel_type = "kEpsilonViscosity";
  auto params = getFactory().getValidParams(kernel_type);
  assignBlocks(params, _blocks);

  params.set<AuxVariableName>("variable") = _turbulent_viscosity_name;
  params.set<MooseFunctorName>(NS::density) = _flow_equations_physics->densityName();
  params.set<MooseFunctorName>(NS::mu) = _flow_equations_physics->dynamicViscosityName();
  params.set<MooseFunctorName>(NS::TKE) = _tke_name;
  params.set<MooseFunctorName>(NS::TKED) = _tked_name;
  params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
  params.set<MooseEnum>("wall_treatment") = _wall_treatment_eps;
  for (const auto d : make_range(dimension()))
    params.set<MooseFunctorName>(u_names[d]) = _velocity_names[d];

  params.applySpecificParameters(parameters(), {"C_mu", "bulk_wall_treatment", "mu_t_ratio_max"});
  params.applySpecificParameters(parameters(), _variant_parameters);
  if (isParamSetByUser("time_scale_Ct"))
    params.set<Real>("Ct") = getParam<Real>("time_scale_Ct");
  if (needsWallDistance())
    params.set<MooseFunctorName>(NS::wall_distance) = _wall_distance_name;
  params.set<ExecFlagEnum>("execute_on") = {EXEC_NONLINEAR};

  getProblem().addAuxKernel(kernel_type, name() + "_viscosity_aux", params);
}

void
AdvancedKEpsilonPhysics::addKEpsilonSink()
{
  const std::string u_names[3] = {"u", "v", "w"};
  const bool use_buoyancy = getParam<bool>("use_buoyancy");

  // Common setup of the two source and sink kernels
  auto set_common_parameters = [&](InputParameters & params)
  {
    assignBlocks(params, _blocks);
    params.set<MooseFunctorName>(NS::density) = _flow_equations_physics->densityName();
    params.set<MooseFunctorName>(NS::mu) = _flow_equations_physics->dynamicViscosityName();
    params.set<MooseFunctorName>(NS::mu_t) = _turbulent_viscosity_name;
    params.set<std::vector<BoundaryName>>("walls") = _turbulence_walls;
    params.set<MooseEnum>("wall_treatment") = _wall_treatment_eps;
    params.set<Real>("C_mu") = getParam<Real>("C_mu");
    params.set<Real>("C_pl") = getParam<Real>("C_pl");
    for (const auto d : make_range(dimension()))
      params.set<MooseFunctorName>(u_names[d]) = _velocity_names[d];

    params.applySpecificParameters(parameters(), _variant_parameters);

    // The buoyancy production term needs the temperature field and the gravity vector of the
    // coupled Physics. The kernels error out if these are specified without buoyancy production
    if (use_buoyancy)
    {
      if (!_has_energy_equation)
        paramError("use_buoyancy",
                   "Buoyancy production of turbulence requires a fluid heat transfer Physics to "
                   "provide the temperature field");
      params.set<MooseFunctorName>("temperature") =
          _fluid_energy_physics->getFluidTemperatureName();
      params.set<MooseFunctorName>("beta") =
          _flow_equations_physics->getParam<MooseFunctorName>("thermal_expansion");
      params.set<RealVectorValue>("gravity") = _flow_equations_physics->gravityVector();
      if (isParamSetByUser(NS::turbulent_Prandtl))
        params.set<Real>(NS::turbulent_Prandtl) = getConstantCoefficient(NS::turbulent_Prandtl);
    }
  };

  {
    const std::string kernel_type = "kEpsilonTKESourceSink";
    InputParameters params = getFactory().getValidParams(kernel_type);
    params.set<LinearVariableName>("variable") = _tke_name;
    params.set<MooseFunctorName>(NS::TKED) = _tked_name;
    set_common_parameters(params);
    getProblem().addLinearFVKernel(kernel_type, prefix() + "tke_source_sink", params);
  }

  {
    const std::string kernel_type = "kEpsilonTKEDSourceSink";
    InputParameters params = getFactory().getValidParams(kernel_type);
    params.set<LinearVariableName>("variable") = _tked_name;
    params.set<MooseFunctorName>(NS::TKE) = _tke_name;
    set_common_parameters(params);
    params.set<Real>("C1_eps") = getConstantCoefficient("C1_eps");
    params.set<Real>("C2_eps") = getConstantCoefficient("C2_eps");
    if (isParamSetByUser("yap_Ct"))
      params.set<Real>("Ct") = getParam<Real>("yap_Ct");
    if (needsWallDistance())
      params.set<MooseFunctorName>(NS::wall_distance) = _wall_distance_name;
    getProblem().addLinearFVKernel(kernel_type, prefix() + "tked_source_sink", params);
  }
}
