//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MoltenSaltRadiolysisAction.h"
#include "FEProblemBase.h"
#include "Factory.h"

#include <algorithm>
#include <sstream>

registerMooseAction("OpenPronghornApp", MoltenSaltRadiolysisAction, "add_variable");
registerMooseAction("OpenPronghornApp", MoltenSaltRadiolysisAction, "add_ic");
registerMooseAction("OpenPronghornApp", MoltenSaltRadiolysisAction, "add_linear_fv_kernel");

namespace
{
// Format a Real with full precision so it can be used as a constant functor name.
std::string
numberToFunctor(Real value)
{
  std::ostringstream oss;
  oss.precision(17);
  oss << value;
  return oss.str();
}
}

InputParameters
MoltenSaltRadiolysisAction::validParams()
{
  InputParameters params = Action::validParams();

  params.addClassDescription("Creates the variables and kernels that integrate a molten salt "
                             "radiolysis problem from the built-in chemistry database.");

  params.addParam<DataFileName>(
      "database",
      "msr_database.json",
      "The JSON chemistry database file (species, reactions, G values, Henry coefficients). "
      "Defaults to the bundled data/msr_database.json; supply your own to change the network.");

  MooseEnum salt_type("chloride fluoride");
  params.addRequiredParam<MooseEnum>("salt_type", salt_type, "The salt kernel to simulate.");

  MooseEnum integration_method("linear epi3v", "linear");
  params.addParam<MooseEnum>(
      "integration_method",
      integration_method,
      "Chemistry integration representation. 'linear' preserves the existing segregated "
      "LinearFV implementation. 'epi3v' uses the coupled state representation required by "
      "MSREPI3VTimeIntegrator.");

  MultiMooseEnum metals("Zn Cr U");
  params.addParam<MultiMooseEnum>("metals", metals, "Metal redox species to include.");

  MultiMooseEnum gases("Cl2 F2");
  params.addParam<MultiMooseEnum>(
      "gas_species", gases, "Diatomic gases tracked with a well-mixed headspace.");

  params.addRequiredParam<MooseFunctorName>("temperature", "The temperature functor [K].");

  params.addParam<MooseFunctorName>(
      "dose_rate",
      "0",
      "Volumetric dose rate functor [J/m^3/s]. May be a constant or a spatially/temporally varying "
      "functor (function, variable, etc.) to represent a non-uniform radiation field.");

  MooseEnum radiation("gamma", "gamma");
  params.addParam<MooseEnum>("radiation", radiation, "The radiation type for the G values.");

  params.addParam<MooseFunctorName>(
      "diffusivity", "Molecular diffusivity for the dissolved species [m^2/s]. Diffusion is added "
                     "only when this is supplied.");

  params.addParam<UserObjectName>(
      "rhie_chow_user_object",
      "The Rhie-Chow user object providing the mass flux. When supplied, the dissolved species are "
      "advected by the segregated flow solution (Navier-Stokes coupling).");

  params.addParam<bool>(
      "use_generic_rhie_chow_advection",
      false,
      "Use OpenPronghorn's MSRLinearFVScalarAdvectionGenericRC instead of MOOSE "
      "LinearFVScalarAdvection. This allows linear FV species systems to consume any "
      "RhieChowFaceFluxProvider, including INSFVRhieChowInterpolator. Default false preserves "
      "legacy behavior.");

  params.addParam<RealVectorValue>(
      "velocity",
      "A prescribed constant advection velocity for the dissolved species. Used only when "
      "'rhie_chow_user_object' is not supplied.");

  MooseEnum interp("average upwind", "upwind");
  params.addParam<MooseEnum>(
      "advected_interp_method", interp, "The advected interpolation method for the species.");

  params.addParam<Real>("kLa", 0.0, "Overall gas-liquid mass-transfer coefficient kLa [1/s].");

  params.addParam<Real>(
      "volume_ratio", 1.0, "Ratio of liquid to headspace volume V_liq / V_gas for gas exchange.");

  params.addParam<RealVectorValue>(
      "gas_buoyancy_velocity",
      "Buoyant rise (slip) velocity of the gas phase relative to the liquid [m/s], typically "
      "directed against gravity. When supplied, the gas-phase species are advected by the flow plus "
      "this slip velocity, capturing buoyant transport of the evolved gas.");

  params.addParam<MooseFunctorName>(
      "gas_dispersivity",
      "Dispersion coefficient for the gas-phase species [m^2/s] (bubble/turbulent dispersion). When "
      "supplied, the gas-phase species are dispersed.");

  params.addParam<std::vector<std::string>>(
      "initial_condition_species", {}, "Species names with a specified initial concentration.");

  params.addParam<std::vector<Real>>("initial_condition_values",
                                     {},
                                     "Initial concentrations [mol/m^3] matching "
                                     "'initial_condition_species' (default 0 for all others).");

  params.addParam<std::vector<std::string>>(
      "g_value_species", {}, "Species names whose default G value is overridden.");

  params.addParam<std::vector<Real>>(
      "g_value_overrides", {}, "G values [molecules/100eV] matching 'g_value_species'.");

  params.addParam<std::vector<SubdomainName>>(
      "block", {}, "The blocks the radiolysis problem is defined on (default: whole mesh).");

  params.addParam<bool>(
      "time_derivative",
      true,
      "Add a time-derivative kernel to every species (transient chemistry). Set to false for a "
      "steady-state dose-driven solve with a Steady or SIMPLE executioner. In epi3v "
      "mode, time advancement will instead be owned by MSREPI3VTimeIntegrator.");

  params.addParam<bool>(
      "transport_only",
      false,
      "In linear integration mode, generate only species variables plus time derivative, "
      "advection, and diffusion kernels. Skip radiolytic source, reaction, and gas-exchange "
      "kernels. This is intended for operator-split transport around an external chemistry "
      "integrator such as EPI3V.");

  params.addParam<bool>("verbose", false, "Print the generated species and solver system names.");

  return params;
}

MoltenSaltRadiolysisAction::MoltenSaltRadiolysisAction(const InputParameters & params)
  : Action(params),
    _salt(getParam<MooseEnum>("salt_type")),
    _integration_method(getParam<MooseEnum>("integration_method")),
    _temperature(getParam<MooseFunctorName>("temperature"))
{
  buildPlan();
}

void
MoltenSaltRadiolysisAction::buildPlan()
{
  // Load the chemistry database from the (input-specified) JSON file.
  MSR::MoltenSaltRadiolysisDatabase db(getParam<DataFileName>("database"));

  // Initial conditions and G value overrides as lookup maps
  const auto & ic_species = getParam<std::vector<std::string>>("initial_condition_species");
  const auto & ic_values = getParam<std::vector<Real>>("initial_condition_values");

  if (ic_species.size() != ic_values.size())
    paramError("initial_condition_values",
               "'initial_condition_species' and 'initial_condition_values' must have equal length.");

  std::map<std::string, Real> ic;
  for (const auto i : index_range(ic_species))
    ic[ic_species[i]] = ic_values[i];

  auto g_values = db.defaultGValues(getParam<MooseEnum>("radiation"), _salt);

  const auto & g_species = getParam<std::vector<std::string>>("g_value_species");
  const auto & g_overrides = getParam<std::vector<Real>>("g_value_overrides");

  if (g_species.size() != g_overrides.size())
    paramError("g_value_overrides",
               "'g_value_species' and 'g_value_overrides' must have equal length.");

  for (const auto i : index_range(g_species))
    g_values[g_species[i]] = g_overrides[i];

  // Assemble the species list: core species, then metal species, then gas-phase species.
  std::vector<std::string> names = db.coreSpecies(_salt);

  for (const auto & metal : getParam<MultiMooseEnum>("metals"))
    for (const auto & sp : db.metalSpecies(metal.name()))
      names.push_back(sp);

  for (const auto & gas : getParam<MultiMooseEnum>("gas_species"))
  {
    GasPair pair;
    pair.liquid = db.gasLiquidSpecies(gas.name());
    pair.vapor = db.gasPhaseSpecies(gas.name());
    pair.kH = db.henryCoefficient(gas.name());
    _gas_pairs.push_back(pair);
    names.push_back(pair.vapor);
  }

  for (const auto & name : names)
  {
    // Gas-phase species are named with a "_gas" suffix; they live in the headspace and only
    // exchange with the melt (no advection, diffusion or radiolytic source).
    const bool is_gas = name.size() > 4 && name.compare(name.size() - 4, 4, "_gas") == 0;

    SpeciesPlan plan;
    plan.name = name;
    plan.is_gas = is_gas;
    plan.initial_condition = ic.count(name) ? ic.at(name) : 0.0;

    const auto g_it = g_values.find(name);
    plan.g_value = (g_it != g_values.end()) ? g_it->second : 0.0;

    _species.push_back(plan);
    _is_gas[name] = is_gas;
  }

  // Assemble the reaction network: core reactions plus the templated metal reactions.
  _reactions = db.coreReactions(_salt);

  for (const auto & metal : getParam<MultiMooseEnum>("metals"))
    for (const auto & rxn : db.metalReactions(metal.name(), _salt))
      _reactions.push_back(rxn);
}

void
MoltenSaltRadiolysisAction::maybeAssignBlocks(InputParameters & params) const
{
  const auto & blocks = getParam<std::vector<SubdomainName>>("block");
  if (!blocks.empty())
    params.set<std::vector<SubdomainName>>("block") = blocks;
}

void
MoltenSaltRadiolysisAction::act()
{
  if (_current_task == "add_variable")
    addVariables();
  else if (_current_task == "add_ic")
    addInitialConditions();
  else if (_current_task == "add_linear_fv_kernel")
    addKernels();
}

void
MoltenSaltRadiolysisAction::addVariables()
{
  if (getParam<bool>("verbose"))
  {
    if (useEPI3V())
    {
      _console << "[MoltenSaltRadiolysis] creating " << _species.size()
               << " species as coupled nonlinear FV variables for EPI3V:\n";
      for (const auto & sp : _species)
        _console << "  " << sp.name << "\n";
    }
    else
    {
      _console << "[MoltenSaltRadiolysis] creating " << _species.size()
               << " species; solver systems (list these in [Problem] linear_sys_names and the "
                  "executioner):\n";
      for (const auto & sp : _species)
        _console << "  " << sp.name << " -> " << systemName(sp.name) << "\n";
    }

    _console << std::flush;
  }

  // The species are finite-volume variables, so the problem needs the finite-volume
  // infrastructure (face information, matrix sparsity); AddVariableAction does this for variables
  // declared in [Variables], so we must do it here when creating them programmatically.
  _problem->needFV();

  for (const auto & sp : _species)
  {
    if (useEPI3V())
    {
      // EPI3V species share the default nonlinear solver system.
      auto params = _factory.getValidParams("MooseVariableFVReal");
      maybeAssignBlocks(params);
      _problem->addVariable("MooseVariableFVReal", sp.name, params);
    }
    else
    {
      auto params = _factory.getValidParams("MooseLinearVariableFVReal");
      params.set<SolverSystemName>("solver_sys") = systemName(sp.name);
      maybeAssignBlocks(params);
      _problem->addVariable("MooseLinearVariableFVReal", sp.name, params);
    }
  }
}

void
MoltenSaltRadiolysisAction::addInitialConditions()
{
  for (const auto & sp : _species)
  {
    auto params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = sp.name;
    params.set<Real>("value") = sp.initial_condition;
    maybeAssignBlocks(params);
    _problem->addInitialCondition("ConstantIC", sp.name + "_ic", params);
  }
}

void
MoltenSaltRadiolysisAction::addKernels()
{
  if (useEPI3V())
  {
    // EPI3V owns chemistry advancement for the coupled nonlinear state.
    if (getParam<bool>("verbose"))
      _console << "[MoltenSaltRadiolysis] epi3v mode: skipping legacy LinearFV chemistry kernels.\n"
               << std::flush;
    return;
  }

  const bool has_rhie_chow = isParamValid("rhie_chow_user_object");
  const bool has_velocity = isParamValid("velocity");
  const bool has_diffusion = isParamValid("diffusivity");
  const bool has_buoyancy = isParamValid("gas_buoyancy_velocity");
  const bool has_gas_dispersion = isParamValid("gas_dispersivity");
  const bool add_time = getParam<bool>("time_derivative");
  const bool transport_only = getParam<bool>("transport_only");

  // Per-species transport and source kernels
  for (const auto & sp : _species)
  {
    // Time derivative for every species (transient chemistry); omitted for steady-state solves.
    if (add_time)
    {
      auto params = _factory.getValidParams("LinearFVTimeDerivative");
      params.set<LinearVariableName>("variable") = sp.name;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVTimeDerivative", "msr_" + sp.name + "_time", params);
    }

    // Gas-phase species: buoyant advection (flow plus rise/slip velocity) and dispersion. They are
    // produced from the melt by the gas-exchange kernels (added below), not by direct radiolysis.
    if (sp.is_gas)
    {
      addGasTransport(sp.name, has_rhie_chow, has_velocity, has_buoyancy, has_gas_dispersion);
      continue;
    }

    if (has_rhie_chow)
    {
      // Advection by the supplied Rhie-Chow face-flux provider.
      // Legacy default uses MOOSE LinearFVScalarAdvection (RhieChowMassFlux-specific).
      // The opt-in generic path accepts any RhieChowFaceFluxProvider, including a reconstructed
      // INSFVRhieChowInterpolator transferred from the steady SIMPLE parent.
      const bool use_generic_rc = getParam<bool>("use_generic_rhie_chow_advection");
      const std::string kernel_type =
          use_generic_rc ? "MSRLinearFVScalarAdvectionGenericRC" : "LinearFVScalarAdvection";

      auto params = _factory.getValidParams(kernel_type);
      params.set<LinearVariableName>("variable") = sp.name;
      params.set<UserObjectName>("rhie_chow_user_object") =
          getParam<UserObjectName>("rhie_chow_user_object");
      params.set<MooseEnum>("advected_interp_method") =
          getParam<MooseEnum>("advected_interp_method");
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel(kernel_type, "msr_" + sp.name + "_advection", params);
    }
    else if (has_velocity)
    {
      // Advection by a prescribed constant velocity field
      auto params = _factory.getValidParams("LinearFVAdvection");
      params.set<LinearVariableName>("variable") = sp.name;
      params.set<RealVectorValue>("velocity") = getParam<RealVectorValue>("velocity");
      params.set<MooseEnum>("advected_interp_method") =
          getParam<MooseEnum>("advected_interp_method");
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVAdvection", "msr_" + sp.name + "_advection", params);
    }

    if (has_diffusion)
    {
      auto params = _factory.getValidParams("LinearFVDiffusion");
      params.set<LinearVariableName>("variable") = sp.name;
      params.set<MooseFunctorName>("diffusion_coeff") = getParam<MooseFunctorName>("diffusivity");
      params.set<bool>("use_nonorthogonal_correction") = false;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVDiffusion", "msr_" + sp.name + "_diffusion", params);
    }

    if (!transport_only && sp.g_value != 0.0)
    {
      // Radiolytic source S_i = G_i * dose_rate / (100 eV * N_A).
      auto params = _factory.getValidParams("LinearFVSource");
      params.set<LinearVariableName>("variable") = sp.name;
      params.set<MooseFunctorName>("source_density") = getParam<MooseFunctorName>("dose_rate");
      params.set<MooseFunctorName>("scaling_factor") =
          numberToFunctor(MSR::gToSource(sp.g_value, 1.0));
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel("LinearFVSource", "msr_" + sp.name + "_radiolysis", params);
    }
  }

  // Operator-split transport is assembled here; EPI3V owns chemistry and source terms.
  if (transport_only)
    return;

  // Reaction kernels
  for (const auto i : index_range(_reactions))
    addReactionKernels(_reactions[i], i);

  // Gas-liquid exchange kernels
  const Real kLa = getParam<Real>("kLa");
  const Real volume_ratio = getParam<Real>("volume_ratio");

  for (const auto & gas : _gas_pairs)
  {
    {
      auto params = _factory.getValidParams("MSRGasExchange");
      params.set<LinearVariableName>("variable") = gas.liquid;
      params.set<MooseEnum>("mode") = "liquid";
      params.set<MooseFunctorName>("partner") = gas.vapor;
      params.set<MooseFunctorName>("temperature") = _temperature;
      params.set<Real>("kLa") = kLa;
      params.set<Real>("kH") = gas.kH;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel(
          "MSRGasExchange", "msr_" + gas.liquid + "_exchange", params);
    }

    {
      auto params = _factory.getValidParams("MSRGasExchange");
      params.set<LinearVariableName>("variable") = gas.vapor;
      params.set<MooseEnum>("mode") = "gas";
      params.set<MooseFunctorName>("partner") = gas.liquid;
      params.set<MooseFunctorName>("temperature") = _temperature;
      params.set<Real>("kLa") = kLa;
      params.set<Real>("kH") = gas.kH;
      params.set<Real>("volume_ratio") = volume_ratio;
      maybeAssignBlocks(params);
      _problem->addLinearFVKernel(
          "MSRGasExchange", "msr_" + gas.vapor + "_exchange", params);
    }
  }
}

void
MoltenSaltRadiolysisAction::addGasTransport(const std::string & name,
                                            bool has_rhie_chow,
                                            bool has_velocity,
                                            bool has_buoyancy,
                                            bool has_gas_dispersion)
{
  const RealVectorValue buoyancy =
      has_buoyancy ? getParam<RealVectorValue>("gas_buoyancy_velocity") : RealVectorValue(0, 0, 0);

  // Advection: by the solved flow plus a buoyant slip velocity (Rhie-Chow), or by a prescribed
  // and/or buoyant constant velocity.
  if (has_rhie_chow)
  {
    const bool use_generic_rc = getParam<bool>("use_generic_rhie_chow_advection");
    const std::string kernel_type =
        use_generic_rc ? "MSRLinearFVScalarAdvectionGenericRC" : "LinearFVScalarAdvection";

    auto params = _factory.getValidParams(kernel_type);
    params.set<LinearVariableName>("variable") = name;
    params.set<UserObjectName>("rhie_chow_user_object") =
        getParam<UserObjectName>("rhie_chow_user_object");
    params.set<MooseEnum>("advected_interp_method") =
        getParam<MooseEnum>("advected_interp_method");

    if (has_buoyancy)
    {
      params.set<MooseFunctorName>("u_slip") = numberToFunctor(buoyancy(0));
      params.set<MooseFunctorName>("v_slip") = numberToFunctor(buoyancy(1));
      params.set<MooseFunctorName>("w_slip") = numberToFunctor(buoyancy(2));
    }

    maybeAssignBlocks(params);
    _problem->addLinearFVKernel(kernel_type, "msr_" + name + "_advection", params);
  }
  else if (has_velocity || has_buoyancy)
  {
    const RealVectorValue flow =
        has_velocity ? getParam<RealVectorValue>("velocity") : RealVectorValue(0, 0, 0);

    auto params = _factory.getValidParams("LinearFVAdvection");
    params.set<LinearVariableName>("variable") = name;
    params.set<RealVectorValue>("velocity") = flow + buoyancy;
    params.set<MooseEnum>("advected_interp_method") =
        getParam<MooseEnum>("advected_interp_method");
    maybeAssignBlocks(params);
    _problem->addLinearFVKernel("LinearFVAdvection", "msr_" + name + "_advection", params);
  }

  // Dispersion of the gas phase (bubble / turbulent dispersion).
  if (has_gas_dispersion)
  {
    auto params = _factory.getValidParams("LinearFVDiffusion");
    params.set<LinearVariableName>("variable") = name;
    params.set<MooseFunctorName>("diffusion_coeff") =
        getParam<MooseFunctorName>("gas_dispersivity");
    params.set<bool>("use_nonorthogonal_correction") = false;
    maybeAssignBlocks(params);
    _problem->addLinearFVKernel("LinearFVDiffusion", "msr_" + name + "_dispersion", params);
  }
}

void
MoltenSaltRadiolysisAction::addReactionKernels(const MSR::ReactionData & reaction,
                                               unsigned int index)
{
  const std::string tag = "msr_rxn" + std::to_string(index) + "_";

  // Helper that copies the Arrhenius parameters and temperature onto a kernel's parameters.
  auto set_kinetics = [&](InputParameters & params)
  {
    params.set<MooseFunctorName>("temperature") = _temperature;
    params.set<Real>("A") = reaction.A;
    params.set<Real>("k_ref") = reaction.k_ref;
    params.set<Real>("T_ref") = reaction.T_ref;
    params.set<Real>("Ea") = reaction.Ea;
  };

  // Consumption: one kernel per reactant term (implicit sink). Partners are the other reactants.
  for (const auto & self : reaction.reactants)
  {
    if (!_is_gas.count(self.species))
      continue; // species not tracked in this problem

    auto params = _factory.getValidParams("MSRReaction");
    params.set<LinearVariableName>("variable") = self.species;
    params.set<MooseEnum>("mode") = "consumption";
    params.set<Real>("stoichiometric_coefficient") = self.coeff;
    params.set<Real>("self_order") = self.order;

    unsigned int n_partner = 0;
    for (const auto & other : reaction.reactants)
    {
      if (&other == &self)
        continue;

      ++n_partner;
      params.set<MooseFunctorName>("partner" + std::to_string(n_partner)) = other.species;
      params.set<Real>("partner" + std::to_string(n_partner) + "_order") = other.order;
    }

    set_kinetics(params);
    maybeAssignBlocks(params);
    _problem->addLinearFVKernel("MSRReaction", tag + self.species + "_consume", params);
  }

  // Production: one kernel per product term (explicit source). Partners are all reactants.
  for (const auto & self : reaction.products)
  {
    if (!_is_gas.count(self.species))
      continue; // species not tracked in this problem

    auto params = _factory.getValidParams("MSRReaction");
    params.set<LinearVariableName>("variable") = self.species;
    params.set<MooseEnum>("mode") = "production";
    params.set<Real>("stoichiometric_coefficient") = self.coeff;

    unsigned int n_partner = 0;
    for (const auto & reactant : reaction.reactants)
    {
      ++n_partner;
      params.set<MooseFunctorName>("partner" + std::to_string(n_partner)) = reactant.species;
      params.set<Real>("partner" + std::to_string(n_partner) + "_order") = reactant.order;
    }

    set_kinetics(params);
    maybeAssignBlocks(params);
    _problem->addLinearFVKernel("MSRReaction", tag + self.species + "_produce", params);
  }
}
