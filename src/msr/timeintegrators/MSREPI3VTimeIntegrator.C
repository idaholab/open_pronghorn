#include "MSREPI3VTimeIntegrator.h"
#include "MSREPI3V.h"

#include "FEProblem.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <atomic>
#include <cstdlib>
#include <thread>
#include <utility>

registerMooseObject("OpenPronghornApp", MSREPI3VTimeIntegrator);

InputParameters
MSREPI3VTimeIntegrator::validParams()
{
  InputParameters params = TimeIntegrator::validParams();

  params.addClassDescription(
      "EPI3V time integrator for cell-local molten-salt radiolysis chemistry on an FV mesh.");

  params.addParam<DataFileName>(
      "database",
      "msr_database.json",
      "JSON chemistry database used to construct the direct MSR chemistry RHS/Jacobian.");

  MooseEnum salt_type("chloride fluoride");
  params.addRequiredParam<MooseEnum>("salt_type", salt_type, "The salt chemistry network.");

  MultiMooseEnum metals("Zn Cr U");
  params.addParam<MultiMooseEnum>("metals", metals, "Metal redox species to include.");

  params.addRequiredParam<std::vector<std::string>>(
      "species",
      "Ordered chemistry species. These names must correspond to nonlinear FV variables.");

  params.addRequiredParam<Real>(
      "temperature",
      "Fallback constant chemistry temperature [K]. Used when temperature_variable is not supplied.");

  params.addParam<std::string>(
      "temperature_variable",
      "",
      "Optional cell-local nonlinear FV temperature variable [K]. When supplied, each chemistry "
      "cell uses this variable instead of the constant temperature parameter.");

  params.addParam<Real>(
      "dose_rate",
      0.0,
      "Fallback constant volumetric dose rate [J/m^3/s]. Used when dose_rate_variable is not supplied.");

  params.addParam<std::string>(
      "dose_rate_variable",
      "",
      "Optional cell-local nonlinear FV volumetric dose-rate variable [J/m^3/s]. When supplied, "
      "each chemistry cell uses this variable instead of the constant dose_rate parameter.");

  MooseEnum radiation("gamma", "gamma");
  params.addParam<MooseEnum>("radiation", radiation, "The radiation type for the G values.");

  params.addParam<std::vector<std::string>>(
      "g_value_species", {}, "Species names whose database G value is overridden.");

  params.addParam<std::vector<Real>>(
      "g_value_overrides",
      {},
      "G values [molecules/100 eV] corresponding to g_value_species.");

  params.addParam<Real>("rtol", 1.0e-6, "Relative tolerance for embedded EPI3V WRMS estimate.");
  params.addParam<Real>("atol", 1.0e-12, "Absolute tolerance for embedded EPI3V WRMS estimate.");

  return params;
}

MSREPI3VTimeIntegrator::MSREPI3VTimeIntegrator(const InputParameters & parameters)
  : TimeIntegrator(parameters),
    _temperature(getParam<Real>("temperature")),
    _temperature_variable(getParam<std::string>("temperature_variable")),
    _dose_rate(getParam<Real>("dose_rate")),
    _dose_rate_variable(getParam<std::string>("dose_rate_variable")),
    _rtol(getParam<Real>("rtol")),
    _atol(getParam<Real>("atol")),
    _last_error_norm(0.0)
{
  if (_temperature <= 0.0)
    paramError("temperature", "Fallback temperature must be positive.");
  if (_dose_rate < 0.0)
    paramError("dose_rate", "Fallback dose_rate must be nonnegative.");
  if (_rtol <= 0.0)
    paramError("rtol", "rtol must be positive.");
  if (_atol <= 0.0)
    paramError("atol", "atol must be positive.");

  _species = getParam<std::vector<std::string>>("species");

  MSR::MoltenSaltRadiolysisDatabase db(getParam<DataFileName>("database"));
  const auto salt_type = getParam<MooseEnum>("salt_type");

  _reactions = db.coreReactions(salt_type);
  for (const auto & metal : getParam<MultiMooseEnum>("metals"))
    for (const auto & reaction : db.metalReactions(metal.name(), salt_type))
      _reactions.push_back(reaction);

  // Apply user G-value overrides to the selected database values.
  auto g_values = db.defaultGValues(getParam<MooseEnum>("radiation"), salt_type);
  const auto & g_species = getParam<std::vector<std::string>>("g_value_species");
  const auto & g_overrides = getParam<std::vector<Real>>("g_value_overrides");

  if (g_species.size() != g_overrides.size())
    paramError("g_value_overrides",
               "'g_value_species' and 'g_value_overrides' must have equal length.");

  for (std::size_t i = 0; i < g_species.size(); ++i)
    g_values[g_species[i]] = g_overrides[i];

  // Apply direct radiolytic sources cell-by-cell in solve().
  _source_per_unit_dose.assign(_species.size(), 0.0);

  for (std::size_t i = 0; i < _species.size(); ++i)
  {
    const auto it = g_values.find(_species[i]);
    if (it != g_values.end())
      _source_per_unit_dose[i] = MSR::gToSource(it->second, 1.0);
  }

  _chemistry = std::make_unique<MSR::MSRChemistrySystem>(_species, _reactions);
}

void
MSREPI3VTimeIntegrator::solve()
{

  auto & mesh = _fe_problem.mesh().getMesh();
  auto & nonlinear_system = *_nonlinear_implicit_system;
  auto & solution = *nonlinear_system.solution;
  const auto & dof_map = nonlinear_system.get_dof_map();

  std::vector<unsigned int> variable_numbers;
  variable_numbers.reserve(_species.size());

  for (const auto & species : _species)
  {
    if (!nonlinear_system.has_variable(species))
      mooseError("MSREPI3VTimeIntegrator: solver system does not contain species variable '",
                 species,
                 "'.");

    variable_numbers.push_back(nonlinear_system.variable_number(species));
  }

  const bool use_local_temperature = !_temperature_variable.empty();
  unsigned int temperature_variable_number = 0;

  if (use_local_temperature)
  {
    if (!nonlinear_system.has_variable(_temperature_variable))
      mooseError("MSREPI3VTimeIntegrator: temperature_variable '",
                 _temperature_variable,
                 "' is not present in the same nonlinear solver system as the chemistry species.");

    temperature_variable_number =
        nonlinear_system.variable_number(_temperature_variable);
  }

  const bool use_local_dose_rate = !_dose_rate_variable.empty();
  unsigned int dose_rate_variable_number = 0;

  if (use_local_dose_rate)
  {
    if (!nonlinear_system.has_variable(_dose_rate_variable))
      mooseError("MSREPI3VTimeIntegrator: dose_rate_variable '",
                 _dose_rate_variable,
                 "' is not present in the same nonlinear solver system as the chemistry species.");

    dose_rate_variable_number =
        nonlinear_system.variable_number(_dose_rate_variable);
  }

  struct PendingCellUpdate
  {
    std::vector<dof_id_type> dofs;
    std::vector<Real> values;
  };

  struct CellWork
  {
    std::vector<dof_id_type> dofs;
    std::vector<Real> y_start;
    Real temperature = 0.0;
    Real dose_rate = 0.0;
  };

  struct WorkerStats
  {
    unsigned int rejected_steps = 0;
    unsigned int accepted_steps = 0;
  };

  std::vector<CellWork> cell_work;

  unsigned int local_rejected_steps = 0;
  unsigned int local_accepted_steps = 0;

  Real local_temperature_min = std::numeric_limits<Real>::infinity();
  Real local_temperature_max = -std::numeric_limits<Real>::infinity();
  Real local_dose_rate_min = std::numeric_limits<Real>::infinity();
  Real local_dose_rate_max = -std::numeric_limits<Real>::infinity();

  // Gather mesh and solution data before the parallel chemistry solve.
  for (const auto * elem : mesh.active_local_element_ptr_range())
  {

    std::vector<dof_id_type> dofs(_species.size());
    std::vector<Real> y_start(_species.size(), 0.0);

    for (std::size_t i = 0; i < _species.size(); ++i)
    {
      std::vector<dof_id_type> variable_dofs;
      dof_map.dof_indices(elem, variable_dofs, variable_numbers[i]);

      if (variable_dofs.size() != 1)
        mooseError("MSREPI3VTimeIntegrator expects exactly one FV dof for species '",
                   _species[i],
                   "' on each active element, but found ",
                   variable_dofs.size(),
                   ".");

      dofs[i] = variable_dofs[0];
      y_start[i] = solution(dofs[i]);
    }

    Real cell_temperature = _temperature;

    if (use_local_temperature)
    {
      std::vector<dof_id_type> temperature_dofs;
      dof_map.dof_indices(elem, temperature_dofs, temperature_variable_number);

      if (temperature_dofs.size() != 1)
        mooseError("MSREPI3VTimeIntegrator expects exactly one FV dof for temperature_variable '",
                   _temperature_variable,
                   "' on each active element, but found ",
                   temperature_dofs.size(),
                   ".");

      cell_temperature = solution(temperature_dofs[0]);
    }

    if (!std::isfinite(cell_temperature) || cell_temperature <= 0.0)
      mooseError("MSREPI3VTimeIntegrator found invalid chemistry temperature ",
                 cell_temperature,
                 " K in element ",
                 elem->id(),
                 ".");

    local_temperature_min = std::min(local_temperature_min, cell_temperature);
    local_temperature_max = std::max(local_temperature_max, cell_temperature);

    Real cell_dose_rate = _dose_rate;

    if (use_local_dose_rate)
    {
      std::vector<dof_id_type> dose_rate_dofs;
      dof_map.dof_indices(elem, dose_rate_dofs, dose_rate_variable_number);

      if (dose_rate_dofs.size() != 1)
        mooseError("MSREPI3VTimeIntegrator expects exactly one FV dof for dose_rate_variable '",
                   _dose_rate_variable,
                   "' on each active element, but found ",
                   dose_rate_dofs.size(),
                   ".");

      cell_dose_rate = solution(dose_rate_dofs[0]);
    }

    if (!std::isfinite(cell_dose_rate) || cell_dose_rate < 0.0)
      mooseError("MSREPI3VTimeIntegrator found invalid volumetric dose rate ",
                 cell_dose_rate,
                 " J/m^3/s in element ",
                 elem->id(),
                 ".");

    local_dose_rate_min = std::min(local_dose_rate_min, cell_dose_rate);
    local_dose_rate_max = std::max(local_dose_rate_max, cell_dose_rate);

    cell_work.push_back(
        {std::move(dofs), std::move(y_start), cell_temperature, cell_dose_rate});
  }
  std::vector<PendingCellUpdate> pending_updates(cell_work.size());

  // Use one chemistry worker unless MSR_EPI3V_THREADS is set.
  unsigned int requested_threads = 1;
  if (const char * env_threads = std::getenv("MSR_EPI3V_THREADS"))
  {
    char * end_ptr = nullptr;
    const unsigned long parsed = std::strtoul(env_threads, &end_ptr, 10);

    if (end_ptr == env_threads || *end_ptr != '\0' || parsed == 0)
      mooseError("MSREPI3VTimeIntegrator: MSR_EPI3V_THREADS must be a positive integer.");

    requested_threads = static_cast<unsigned int>(parsed);
  }

  const unsigned int chemistry_threads =
      cell_work.empty()
          ? 1u
          : std::max(1u,
                     std::min(requested_threads,
                              static_cast<unsigned int>(cell_work.size())));

  std::vector<WorkerStats> worker_stats(chemistry_threads);
  std::atomic<std::size_t> next_cell{0};
  const Real negativity_tolerance = 10.0 * _atol;

  // Advance independent cell-local chemistry in parallel.

  const auto worker = [&](const unsigned int worker_id)
  {
    WorkerStats & stats = worker_stats[worker_id];

    // Each worker owns its chemistry system.
    MSR::MSRChemistrySystem chemistry(_species, _reactions);

    while (true)
    {
      const std::size_t cell_index =
          next_cell.fetch_add(1, std::memory_order_relaxed);

      if (cell_index >= cell_work.size())
        break;

      const auto & work = cell_work[cell_index];

      const MSR::MSREPI3V::RHSFunction rhs =
          [this, &chemistry, &work](const std::vector<Real> & state)
          {

            auto value = chemistry.evaluateRHS(state, work.temperature);

            for (std::size_t i = 0; i < value.size(); ++i)
              value[i] += _source_per_unit_dose[i] * work.dose_rate;
            return value;
          };

      const MSR::MSREPI3V::JacobianFunction jacobian =
          [&chemistry, &work](const std::vector<Real> & state)
          {

            auto value = chemistry.evaluateJacobian(state, work.temperature);
            return value;
          };

      const auto adaptive = MSR::MSREPI3V::integrateAdaptive(
          work.y_start,
          0.0,
          _dt,
          std::min(_dt, std::max(1.0e-12, 0.1 * _dt)),
          _rtol,
          _atol,
          rhs,
          jacobian,
          true,
          negativity_tolerance);

      stats.rejected_steps += adaptive.rejected_steps;
      stats.accepted_steps += adaptive.accepted_steps;

      pending_updates[cell_index] = {work.dofs, adaptive.solution};
    }
  };

  if (chemistry_threads == 1)
    worker(0);
  else
  {
    std::vector<std::thread> workers;
    workers.reserve(chemistry_threads);

    for (unsigned int worker_id = 0; worker_id < chemistry_threads; ++worker_id)
      workers.emplace_back(worker, worker_id);

    for (auto & thread : workers)
      thread.join();
  }

  for (const auto & stats : worker_stats)
  {
    local_rejected_steps += stats.rejected_steps;
    local_accepted_steps += stats.accepted_steps;
  }

  unsigned int global_rejected_steps = local_rejected_steps;
  unsigned int global_accepted_steps = local_accepted_steps;
  comm().sum(global_rejected_steps);
  comm().sum(global_accepted_steps);

  Real global_temperature_min = local_temperature_min;
  Real global_temperature_max = local_temperature_max;
  comm().min(global_temperature_min);
  comm().max(global_temperature_max);

  Real global_dose_rate_min = local_dose_rate_min;
  Real global_dose_rate_max = local_dose_rate_max;
  comm().min(global_dose_rate_min);
  comm().max(global_dose_rate_max);

  for (const auto & update : pending_updates)
    for (std::size_t i = 0; i < update.values.size(); ++i)
      solution.set(update.dofs[i], update.values[i]);

  solution.close();
  nonlinear_system.update();

  _last_error_norm = 0.0;

  _console << "EPI3V internal adaptive chemistry: accepted=" << global_accepted_steps
           << ", rejected=" << global_rejected_steps;

  if (use_local_temperature)
    _console << ", temperature_range_K=[" << global_temperature_min << ", "
             << global_temperature_max << "]";

  if (use_local_dose_rate)
    _console << ", dose_rate_range_J_m3_s=[" << global_dose_rate_min << ", "
             << global_dose_rate_max << "]";

  _console << std::endl;


  if (_nl)
    _nl->setSolution(*nonlinear_system.current_local_solution);

  _n_nonlinear_iterations = 0;
  _n_linear_iterations = 0;

  if (nonlinear_system.nonlinear_solver)
    nonlinear_system.nonlinear_solver->converged = true;
}

void
MSREPI3VTimeIntegrator::computeTimeDerivatives()
{
  if (!_sys.solutionUDot())
    mooseError("MSREPI3VTimeIntegrator: solution u_dot storage was not allocated.");

  NumericVector<Number> & u_dot = *_sys.solutionUDot();
  u_dot = *_solution;
  u_dot -= _solution_old;
  u_dot *= 1.0 / _dt;
  u_dot.close();

  computeDuDotDu();
}

void
MSREPI3VTimeIntegrator::computeADTimeDerivatives(ADReal & ad_u_dot,
                                                 const dof_id_type & dof,
                                                 ADReal & /* ad_u_dotdot */) const
{
  ad_u_dot -= _solution_old(dof);
  ad_u_dot *= 1.0 / _dt;
}

void
MSREPI3VTimeIntegrator::postResidual(NumericVector<Number> & residual)
{
  residual.zero();
  residual.close();
}
