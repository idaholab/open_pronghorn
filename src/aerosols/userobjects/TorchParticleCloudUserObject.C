//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#if defined(MOOSE_LIBTORCH_ENABLED) || defined(LIBTORCH_ENABLED)

#include "TorchParticleCloudUserObject.h"

#include <algorithm>
#include <cmath>
#include <random>

registerMooseObject("OpenPronghornApp", TorchParticleCloudUserObject);

const libMesh::VectorValue<Real> TorchParticleCloudUserObject::_zero_vec(0.0, 0.0, 0.0);

InputParameters
TorchParticleCloudUserObject::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += ADFunctorInterface::validParams();

  params.addClassDescription("MPI-capable torch Lagrangian particle cloud with stable drag update and robust migration.");

  // Fluid functors
  params.addRequiredParam<MooseFunctorName>("u", "x-velocity functor");
  params.addParam<MooseFunctorName>("v", "y-velocity functor (required for 2D/3D)");
  params.addParam<MooseFunctorName>("w", "z-velocity functor (required for 3D)");
  params.addRequiredParam<MooseFunctorName>("rho", "fluid density functor");
  params.addRequiredParam<MooseFunctorName>("mu", "fluid dynamic viscosity functor");

  // Particles (GLOBAL)
  params.addParam<unsigned int>("num_particles", 0, "Initial GLOBAL number of computational particles/parcels (0 allowed when using injection).");
  params.addRequiredParam<Real>("particle_density", "Particle density (kg/m^3)");
  params.addRequiredParam<Real>("particle_diameter", "Particle diameter (m)");
  params.addParam<Real>("parcel_weight", 1.0, "Physical particles represented per computational particle.");

  params.addParam<Point>("init_box_min", Point(0,0,0), "Min corner of uniform seeding box");
  params.addParam<Point>("init_box_max", Point(0,0,0), "Max corner of uniform seeding box");
  params.addParam<RealVectorValue>("init_velocity", RealVectorValue(0,0,0), "Initial particle velocity (m/s)");


  // Injection (optional)
  params.addParam<Real>("injection_rate", 0.0,
                      "Physical particle injection rate [#/s] (global). Converted to parcels using parcel_weight.");
  params.addParam<Real>("injection_rate_parcels", 0.0,
                      "Computational parcel injection rate [parcels/s] (global). If >0, overrides injection_rate.");
  params.addParam<Real>("injection_start_time", 0.0, "Injection start time [s]");
  params.addParam<Real>("injection_end_time", -1.0, "Injection end time [s]. <0 => no end.");
  params.addParam<Point>("injection_box_min", Point(0,0,0), "Min corner of uniform injection box");
  params.addParam<Point>("injection_box_max", Point(0,0,0), "Max corner of uniform injection box");
  params.addParam<RealVectorValue>("injection_velocity", RealVectorValue(0,0,0), "Initial velocity of injected parcels (m/s)");

  params.addParam<RealVectorValue>("gravity", RealVectorValue(0,0,-9.81), "Gravity acceleration (m/s^2)");

  params.addParam<unsigned int>("max_substeps", 50, "Max particle substeps per FV step");
  params.addParam<Real>("max_cfl", 0.5, "Particle CFL based on min cell length (<=0 disables)");

  params.addParam<unsigned int>("rng_seed", 1234, "Base RNG seed (rank-unique seed derived from this)");
  params.addParam<bool>("debug", false, "Print global mapped count + max conc + vmax once per time value.");

  // MPI behavior
  params.addParam<bool>("drop_out_of_mesh", true, "Drop particles that cannot be located / routed.");
  params.addParam<bool>("bbox_fallback", true, "If PointLocator fails, route by per-rank bounding boxes (collective allgather).");
  params.addParam<int>("migration_tag_base", 23000, "Base MPI tag block (avoid collisions with other objects).");

  // Debug synchronization / root-serialized migration
  params.addParam<bool>("mpi_barriers", false, "Insert MPI barriers between particle stages (debug aid; can be slow).");
  params.addParam<bool>("migration_via_root", true,
                        "Serialize migration through a root rank (robust debugging baseline; avoids all-to-all deadlocks).");
  params.addParam<unsigned int>("root_rank", 0, "Root rank used for migration_via_root.");

  // Safety clamps
  params.addParam<Real>("mu_min", 1e-10, "Minimum viscosity used in drag computations [Pa*s]");
  params.addParam<Real>("re_max", 1e6, "Maximum Reynolds number used in Cd correlation");
  params.addParam<Real>("cd_max", 2e3, "Maximum Cd used (safety clamp)");
  params.addParam<Real>("k_max", 1e9, "Maximum drag rate (1/s) used in exponential update");
  params.addParam<Real>("v_max", 1e3, "Maximum particle speed magnitude (m/s) safeguard");


  // Wall deposition (optional)
  params.addParam<Real>("wall_deposition_velocity", 0.0,
                      "Effective wall deposition velocity Vd [m/s]. 0 disables wall deposition sink.");
  params.addParam<bool>("include_settling_in_vd", false,
                      "If true, add an approximate Stokes settling velocity to Vd (rough, isotropic).");

  params.set<ExecFlagEnum>("execute_on", true) = {EXEC_TIMESTEP_END};

  return params;
}

TorchParticleCloudUserObject::TorchParticleCloudUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    ADFunctorInterface(this),

    _mesh_dim(_subproblem.mesh().dimension()),

    _num_particles_initial(getParam<unsigned int>("num_particles")),
    _particle_density(getParam<Real>("particle_density")),
    _particle_diameter(getParam<Real>("particle_diameter")),
    _parcel_weight(getParam<Real>("parcel_weight")),

    _init_box_min(getParam<Point>("init_box_min")),
    _init_box_max(getParam<Point>("init_box_max")),
    _init_velocity(getParam<RealVectorValue>("init_velocity")),
    _enable_injection((getParam<Real>("injection_rate_parcels") > 0.0) || (getParam<Real>("injection_rate") > 0.0)),
    _injection_rate_parcels((getParam<Real>("injection_rate_parcels") > 0.0)
                          ? getParam<Real>("injection_rate_parcels")
                          : (getParam<Real>("injection_rate") / std::max(getParam<Real>("parcel_weight"), 1e-30))),
    _injection_start_time(getParam<Real>("injection_start_time")),
    _injection_end_time(getParam<Real>("injection_end_time")),
    _injection_box_min(getParam<Point>("injection_box_min")),
    _injection_box_max(getParam<Point>("injection_box_max")),
    _injection_velocity(getParam<RealVectorValue>("injection_velocity")),

    _wall_deposition_velocity(getParam<Real>("wall_deposition_velocity")),
    _include_settling_in_vd(getParam<bool>("include_settling_in_vd")),


    _gravity(getParam<RealVectorValue>("gravity")),
    _max_substeps(getParam<unsigned int>("max_substeps")),
    _max_cfl(getParam<Real>("max_cfl")),

    _rng_seed(getParam<unsigned int>("rng_seed")),
    _debug(getParam<bool>("debug")),

    _drop_out_of_mesh(getParam<bool>("drop_out_of_mesh")),
    _bbox_fallback(getParam<bool>("bbox_fallback")),
    _migration_tag_base(getParam<int>("migration_tag_base")),

    _mpi_barriers(getParam<bool>("mpi_barriers")),
    _migration_via_root(getParam<bool>("migration_via_root")),
    _root_rank(getParam<unsigned int>("root_rank")),

    _mu_min(getParam<Real>("mu_min")),
    _re_max(getParam<Real>("re_max")),
    _cd_max(getParam<Real>("cd_max")),
    _k_max(getParam<Real>("k_max")),
    _v_max(getParam<Real>("v_max")),

    _u(getFunctor<ADReal>("u")),
    _v(isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _w(isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr),
    _rho(getFunctor<ADReal>("rho")),
    _mu(getFunctor<ADReal>("mu")),

    _device(_app.getLibtorchDevice())
{
  if (_mesh_dim > 1 && !_v)
    mooseError("TorchParticleCloudUserObject: mesh_dim > 1 but no 'v' functor provided.");
  if (_mesh_dim > 2 && !_w)
    mooseError("TorchParticleCloudUserObject: mesh_dim > 2 but no 'w' functor provided.");


  // Require either an initial population or an injection source
  if (_num_particles_initial == 0 && (!_enable_injection || _injection_rate_parcels <= 0.0))
    mooseError("TorchParticleCloudUserObject: either set 'num_particles' > 0 or enable injection via 'injection_rate'/'injection_rate_parcels'.");

  if (_injection_rate_parcels < 0.0)
    mooseError("TorchParticleCloudUserObject: injection_rate/injection_rate_parcels must be >= 0.");

  if (_injection_end_time >= 0.0 && _injection_end_time < _injection_start_time)
    mooseError("TorchParticleCloudUserObject: injection_end_time must be >= injection_start_time (or <0 to disable end time).");

  if (_wall_deposition_velocity < 0.0)
    mooseError("TorchParticleCloudUserObject: wall_deposition_velocity must be >= 0.");

  if (_particle_density <= 0.0 || _particle_diameter <= 0.0 || _parcel_weight <= 0.0)
    mooseError("TorchParticleCloudUserObject: particle_density/particle_diameter/parcel_weight must be > 0.");
  if (_mu_min <= 0.0 || _re_max <= 0.0 || _cd_max <= 0.0 || _k_max <= 0.0 || _v_max <= 0.0)
    mooseError("TorchParticleCloudUserObject: mu_min/re_max/cd_max/k_max/v_max must be > 0.");

  // Rank info
  const libMesh::MeshBase & mesh = _subproblem.mesh().getMesh();
  _rank = mesh.processor_id();
  _n_procs = mesh.n_processors();

  if (_root_rank >= _n_procs)
    mooseError("TorchParticleCloudUserObject: root_rank must be < n_processors.");

  // Local particle count by simple block distribution (global deterministic)
  const unsigned int base = _num_particles_initial / _n_procs;
  const unsigned int rem  = _num_particles_initial % _n_procs;
  _n_local = static_cast<long>(base + ((_rank < rem) ? 1u : 0u));

  // Rank-unique deterministic seed
  const unsigned int prime = 104729u;
  torch::manual_seed(_rng_seed + prime * static_cast<unsigned int>(_rank));
}

void
TorchParticleCloudUserObject::initialize()
{
  // No MPI comm here (avoid deadlocks across MOOSE phases)
  if (!_cache_built)
    buildCellCache();

  if (!_x.defined())
  {
    auto opts = torch::TensorOptions().dtype(torch::kFloat64).device(_device);

    auto lo = torch::tensor({_init_box_min(0), _init_box_min(1), _init_box_min(2)}, opts);
    auto hi = torch::tensor({_init_box_max(0), _init_box_max(1), _init_box_max(2)}, opts);

    const long N = _n_local;

    _x = torch::rand({N, 3}, opts) * (hi - lo) + lo;

    _v_p = torch::zeros({N, 3}, opts);
    if (N > 0)
    {
      _v_p.index_put_({torch::indexing::Slice(), 0}, _init_velocity(0));
      _v_p.index_put_({torch::indexing::Slice(), 1}, _init_velocity(1));
      _v_p.index_put_({torch::indexing::Slice(), 2}, _init_velocity(2));
    }

    _omega = torch::zeros({N, 3}, opts);

    const Real volume = (libMesh::pi / 6.0) * std::pow(_particle_diameter, 3);
    const Real m = _particle_density * volume;
    _mass = torch::full({N, 1}, m, opts);

    _cell_idx  = torch::full({N}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));
    _dest_rank = torch::full({N}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));
    _f_drag    = torch::zeros({N, 3}, opts);
  }
}

void
TorchParticleCloudUserObject::finalize()
{
}

void
TorchParticleCloudUserObject::execute()
{
  if (!_cache_built)
    buildCellCache();

  ++_step_counter;

  auto & comm = _subproblem.mesh().getMesh().comm();
  if (_mpi_barriers && _n_procs > 1)
    comm.barrier();

  ensureProcessorBoundingBoxes();

  if (_mpi_barriers && _n_procs > 1)
    comm.barrier();

  // Clear outputs
  _cell_forces.clear();
  std::fill(_cell_number_conc.begin(), _cell_number_conc.end(), 0.0);
  std::fill(_cell_mass_conc.begin(), _cell_mass_conc.end(), 0.0);

  sampleFluidToCells();

  // Optional rate-based injection (global deterministic count, locally seeded positions)
  injectParticles();

  if (_mpi_barriers && _n_procs > 1)
    comm.barrier();

  // Pre-step migrate so particles reside on owning rank (as best we can)
  mapParticlesToCells();
  migrateParticles();
  mapParticlesToCells();

  if (_mpi_barriers && _n_procs > 1)
    comm.barrier();

  // Advance
  advanceParticlesExplicit();

  if (_mpi_barriers && _n_procs > 1)
    comm.barrier();

  // Post-step migrate and map for coupling
  mapParticlesToCells();
  migrateParticles();
  mapParticlesToCells();

  if (_mpi_barriers && _n_procs > 1)
    comm.barrier();

// Accumulate reaction forces on local owned cells (drag coupling)
  if (_n_local > 0 && _f_drag.defined() && _f_drag.numel() > 0)
  {
    auto f_cpu = _f_drag.to(torch::kCPU).contiguous();
    auto idx_cpu = _cell_idx.to(torch::kCPU).contiguous();

    auto f = f_cpu.accessor<double, 2>();
    auto idx = idx_cpu.accessor<long, 1>();

    for (long i = 0; i < _n_local; ++i)
    {
      const long c = idx[i];
      if (c < 0)
        continue;
      if (static_cast<std::size_t>(c) >= _local_elems.size())
        continue;

      const libMesh::Elem * elem = _local_elems[static_cast<std::size_t>(c)];
      if (!elem)
        continue;

      const dof_id_type eid = elem->id();

      libMesh::VectorValue<Real> rf(-_parcel_weight * static_cast<Real>(f[i][0]),
                                    -_parcel_weight * static_cast<Real>(f[i][1]),
                                    -_parcel_weight * static_cast<Real>(f[i][2]));

      auto it = _cell_forces.find(eid);
      if (it == _cell_forces.end())
        _cell_forces.emplace(eid, rf);
      else
        it->second += rf;
    }
  }

// Optional wall deposition sink: remove a fraction of parcels in wall-adjacent cells
applyWallDeposition();

// Concentrations are computed after deposition so they reflect the "end-of-step" cloud state
computeCellConcentrations();

// Debug reductions MUST be executed by all ranks if enabled
  if (_debug && (!(_last_debug_time == _t)))
  {
    _last_debug_time = _t;

    auto valid = _cell_idx.ge(0);
    const auto n_valid_local = valid.sum().item<int64_t>();

    Real max_n_local = 0.0, max_m_local = 0.0;
    for (unsigned int c = 0; c < _cell_number_conc.size(); ++c)
    {
      max_n_local = std::max(max_n_local, _cell_number_conc[c]);
      max_m_local = std::max(max_m_local, _cell_mass_conc[c]);
    }

    Real vmax_local = 0.0;
    std::size_t bad_local = 0;
    {
      auto v_cpu = _v_p.to(torch::kCPU).contiguous();
      auto vacc = v_cpu.accessor<double, 2>();
      for (long i = 0; i < v_cpu.size(0); ++i)
      {
        const double vx = vacc[i][0], vy = vacc[i][1], vz = vacc[i][2];
        if (!std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz))
          continue;
        const double s = std::sqrt(vx*vx + vy*vy + vz*vz);
        vmax_local = std::max(vmax_local, static_cast<Real>(s));
      }

      auto x_cpu = _x.to(torch::kCPU).contiguous();
      auto xacc = x_cpu.accessor<double, 2>();
      for (long i = 0; i < x_cpu.size(0); ++i)
        for (int d = 0; d < 3; ++d)
          bad_local += !std::isfinite(xacc[i][d]);
    }

    int64_t n_valid_global = n_valid_local;
    comm.sum(n_valid_global);

    int64_t n_total_global = static_cast<int64_t>(_n_local);
    comm.sum(n_total_global);

    unsigned long long dep_step_global = _n_deposited_step;
    unsigned long long dep_total_global = _n_deposited_total;
    comm.sum(dep_step_global);
    comm.sum(dep_total_global);

    Real max_n = max_n_local, max_m = max_m_local, vmax = vmax_local;
    comm.max(max_n);
    comm.max(max_m);
    comm.max(vmax);

    std::size_t bad_global = bad_local;
    comm.sum(bad_global);

    if (_rank == _root_rank)
      mooseInfoRepeated("TorchParticleCloudUserObject(MPI): mapped=", n_valid_global, "/", n_total_global,
                        " | max nconc=", max_n, " [#/m^3], max mconc=", max_m, " [kg/m^3]",
                        " | vmax=", vmax, " [m/s] | nonfinite(x)=", bad_global,
                        " | deposited(step/total)=", dep_step_global, "/", dep_total_global);
  }

  if (_mpi_barriers && _n_procs > 1)
    comm.barrier();
}

const libMesh::VectorValue<Real> &
TorchParticleCloudUserObject::cellForce(const libMesh::Elem * elem) const
{
  if (!elem)
    return _zero_vec;

  const auto it = _cell_forces.find(elem->id());
  if (it == _cell_forces.end())
    return _zero_vec;
  return it->second;
}

libMesh::VectorValue<Real>
TorchParticleCloudUserObject::cellForceDensity(const libMesh::Elem * elem) const
{
  if (!elem)
    return _zero_vec;

  const auto it_map = _cell_forces.find(elem->id());
  if (it_map == _cell_forces.end())
    return _zero_vec;

  const auto it_loc = _elem_id_to_local.find(elem->id());
  if (it_loc == _elem_id_to_local.end())
    return _zero_vec;

  const unsigned int local_idx = it_loc->second;
  if (local_idx >= _cell_volume.size())
    return _zero_vec;

  const Real V = _cell_volume[local_idx];
  if (V <= 0.0)
    return _zero_vec;

  return it_map->second / V;
}

Real
TorchParticleCloudUserObject::cellNumberConcentration(const libMesh::Elem * elem) const
{
  if (!elem)
    return 0.0;
  const auto it = _elem_id_to_local.find(elem->id());
  if (it == _elem_id_to_local.end())
    return 0.0;
  const unsigned int idx = it->second;
  return (idx < _cell_number_conc.size()) ? _cell_number_conc[idx] : 0.0;
}

Real
TorchParticleCloudUserObject::cellMassConcentration(const libMesh::Elem * elem) const
{
  if (!elem)
    return 0.0;
  const auto it = _elem_id_to_local.find(elem->id());
  if (it == _elem_id_to_local.end())
    return 0.0;
  const unsigned int idx = it->second;
  return (idx < _cell_mass_conc.size()) ? _cell_mass_conc[idx] : 0.0;
}

void
TorchParticleCloudUserObject::buildCellCache()
{
  _local_elems.clear();
  _cell_centroid.clear();
  _cell_volume.clear();
  _cell_wall_area.clear();
  _cell_vel.clear();
  _cell_rho.clear();
  _cell_mu.clear();
  _elem_id_to_local.clear();
  _min_cell_h = std::numeric_limits<Real>::max();

  const libMesh::MeshBase & mesh = _subproblem.mesh().getMesh();

  for (auto it = mesh.active_local_elements_begin(); it != mesh.active_local_elements_end(); ++it)
  {
    const libMesh::Elem * elem = *it;
    if (!elem)
      continue;

    const unsigned int idx = static_cast<unsigned int>(_local_elems.size());
    _local_elems.push_back(elem);
    _elem_id_to_local.emplace(elem->id(), idx);

    _cell_centroid.push_back(elem->true_centroid());
    const Real V = elem->volume();
    _cell_volume.push_back(V);

    // Boundary measure for simple wall-deposition models
    Real A_wall = 0.0;
    for (unsigned int s = 0; s < elem->n_sides(); ++s)
    {
      if (elem->neighbor_ptr(s) == nullptr)
      {
        auto side = elem->side_ptr(s);
        if (side)
          A_wall += side->volume(); // area in 3D, edge length in 2D
      }
    }
    _cell_wall_area.push_back(A_wall);

    const Real h = std::pow(std::max(V, 1e-30), 1.0 / std::max(1u, _mesh_dim));
    _min_cell_h = std::min(_min_cell_h, h);

    _cell_vel.emplace_back(RealVectorValue(0,0,0));
    _cell_rho.emplace_back(0.0);
    _cell_mu.emplace_back(0.0);
  }

  _point_locator = _subproblem.mesh().getMesh().sub_point_locator();
  _point_locator->enable_out_of_mesh_mode();

  _cell_number_conc.assign(_local_elems.size(), 0.0);
  _cell_mass_conc.assign(_local_elems.size(), 0.0);

  if (!_cell_volume.empty())
  {
      const long Nc = static_cast<long>(_cell_volume.size());
      auto cpu = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);

      _cell_volume_t = torch::from_blob(_cell_volume.data(), {Nc}, cpu).clone().to(_device);

      // Precompute (A_wall / V_cell) on device for deposition models
      std::vector<double> aov(static_cast<std::size_t>(Nc), 0.0);
      for (long c = 0; c < Nc; ++c)
        aov[static_cast<std::size_t>(c)] = static_cast<double>(_cell_wall_area[static_cast<std::size_t>(c)]) /
                                           (static_cast<double>(_cell_volume[static_cast<std::size_t>(c)]) + 1e-30);

      _cell_aw_over_v_t = torch::from_blob(aov.data(), {Nc}, cpu).clone().to(_device);
  }
  else
  {
      _cell_volume_t = torch::Tensor();
      _cell_aw_over_v_t = torch::Tensor();
  }

  _cache_built = true;
}

void
TorchParticleCloudUserObject::sampleFluidToCells()
{
  const auto state = determineState();

  for (unsigned int i = 0; i < _local_elems.size(); ++i)
  {
    const libMesh::Elem * elem = _local_elems[i];
    const auto elem_arg = makeElemArg(elem);

    const Real ux = MetaPhysicL::raw_value(_u(elem_arg, state));
    Real uy = 0.0, uz = 0.0;
    if (_mesh_dim > 1)
      uy = MetaPhysicL::raw_value((*_v)(elem_arg, state));
    if (_mesh_dim > 2)
      uz = MetaPhysicL::raw_value((*_w)(elem_arg, state));
    _cell_vel[i] = RealVectorValue(ux, uy, uz);

    const Real rho_val = MetaPhysicL::raw_value(_rho(elem_arg, state));
    const Real mu_val  = MetaPhysicL::raw_value(_mu(elem_arg, state));

    _cell_rho[i] = std::max(rho_val, 1e-30);
    _cell_mu[i]  = std::max(mu_val,  _mu_min);
  }
}

void
TorchParticleCloudUserObject::ensureProcessorBoundingBoxes()
{
  if (!_bbox_fallback || _n_procs <= 1)
    return;

  if (_proc_bbox_built)
    return;

  // local bbox over owned elements
  libMesh::Point mn( std::numeric_limits<Real>::max(),
                     std::numeric_limits<Real>::max(),
                     std::numeric_limits<Real>::max());
  libMesh::Point mx(-std::numeric_limits<Real>::max(),
                    -std::numeric_limits<Real>::max(),
                    -std::numeric_limits<Real>::max());

  for (const auto * e : _local_elems)
  {
    if (!e)
      continue;
    const auto bb = e->loose_bounding_box();
    for (unsigned d = 0; d < 3; ++d)
    {
      mn(d) = std::min(mn(d), bb.min()(d));
      mx(d) = std::max(mx(d), bb.max()(d));
    }
  }

  std::vector<Real> local(6);
  local[0]=mn(0); local[1]=mn(1); local[2]=mn(2);
  local[3]=mx(0); local[4]=mx(1); local[5]=mx(2);

  auto & comm = _subproblem.mesh().getMesh().comm();
  std::vector<Real> all = local;
  comm.allgather(all, /*identical_size=*/true);

  if (all.size() != 6 * static_cast<std::size_t>(_n_procs))
    mooseError("TorchParticleCloudUserObject: unexpected allgather size for bbox exchange.");

  _proc_bbox_min.assign(_n_procs, libMesh::Point());
  _proc_bbox_max.assign(_n_procs, libMesh::Point());

  for (unsigned int r = 0; r < _n_procs; ++r)
  {
    const std::size_t base = 6u * r;
    _proc_bbox_min[r] = libMesh::Point(all[base+0], all[base+1], all[base+2]);
    _proc_bbox_max[r] = libMesh::Point(all[base+3], all[base+4], all[base+5]);
  }

  _proc_bbox_built = true;
}

libMesh::processor_id_type
TorchParticleCloudUserObject::bboxRoute(const libMesh::Point & p) const
{
  if (!_proc_bbox_built)
    return static_cast<libMesh::processor_id_type>(-1);

  const Real eps = 1e-12;
  for (unsigned int r = 0; r < _n_procs; ++r)
  {
    const auto & mn = _proc_bbox_min[r];
    const auto & mx = _proc_bbox_max[r];

    if (!(mn(0) <= mx(0) && mn(1) <= mx(1) && mn(2) <= mx(2)))
      continue;

    if (p(0) >= mn(0)-eps && p(0) <= mx(0)+eps &&
        p(1) >= mn(1)-eps && p(1) <= mx(1)+eps &&
        p(2) >= mn(2)-eps && p(2) <= mx(2)+eps)
      return static_cast<libMesh::processor_id_type>(r);
  }
  return static_cast<libMesh::processor_id_type>(-1);
}

void
TorchParticleCloudUserObject::mapParticlesToCells()
{
  if (!_point_locator)
  {
    _point_locator = _subproblem.mesh().getMesh().sub_point_locator();
    _point_locator->enable_out_of_mesh_mode();
  }

  auto x_cpu = _x.to(torch::kCPU).contiguous();
  auto acc = x_cpu.accessor<double, 2>();

  std::vector<long> cell_host(_n_local, -1);
  std::vector<long> dest_host(_n_local, -1);

  for (long i = 0; i < _n_local; ++i)
  {
    const double px = acc[i][0], py = acc[i][1], pz = acc[i][2];
    if (!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz))
      continue;

    const libMesh::Point p(px, py, pz);

    const libMesh::Elem * elem = (*_point_locator)(p);

    if (elem)
    {
      dest_host[i] = static_cast<long>(elem->processor_id());

      const auto it = _elem_id_to_local.find(elem->id());
      cell_host[i] = (it != _elem_id_to_local.end()) ? static_cast<long>(it->second) : -1;
      continue;
    }

    if (_bbox_fallback && _proc_bbox_built)
    {
      const auto r = bboxRoute(p);
      dest_host[i] = static_cast<long>(r);
      cell_host[i] = -1;
    }
  }

  _cell_idx = torch::from_blob(cell_host.data(), {_n_local}, torch::TensorOptions().dtype(torch::kInt64)).clone().to(_device);
  _dest_rank = torch::from_blob(dest_host.data(), {_n_local}, torch::TensorOptions().dtype(torch::kInt64)).clone().to(_device);
}

void
TorchParticleCloudUserObject::migrateParticles()
{
  if (_n_procs <= 1)
    return;

  if (_migration_via_root)
    migrateParticlesViaRoot();
  else
    migrateParticlesViaRoot(); // for now, keep only the robust mode
}

void
TorchParticleCloudUserObject::migrateParticlesViaRoot()
{
  auto & comm = _subproblem.mesh().getMesh().comm();

  const libMesh::Parallel::MessageTag tag_out_count(_migration_tag_base + 0);
  const libMesh::Parallel::MessageTag tag_out_data (_migration_tag_base + 1);
  const libMesh::Parallel::MessageTag tag_in_count (_migration_tag_base + 2);
  const libMesh::Parallel::MessageTag tag_in_data  (_migration_tag_base + 3);

  // Payload stride: dest + x3 + v3 + omega3 + f_drag3 + mass
  constexpr unsigned int STRIDE = 14;

  // Pull tensors to CPU
  auto x_cpu = _x.to(torch::kCPU).contiguous();
  auto v_cpu = _v_p.to(torch::kCPU).contiguous();
  auto o_cpu = _omega.to(torch::kCPU).contiguous();
  auto f_cpu = _f_drag.to(torch::kCPU).contiguous();
  auto m_cpu = _mass.to(torch::kCPU).contiguous();
  auto d_cpu = _dest_rank.to(torch::kCPU).contiguous();

  auto xa = x_cpu.accessor<double, 2>();
  auto va = v_cpu.accessor<double, 2>();
  auto oa = o_cpu.accessor<double, 2>();
  auto fa = f_cpu.accessor<double, 2>();
  auto ma = m_cpu.accessor<double, 2>();
  auto da = d_cpu.accessor<long, 1>();

  // Pack outgoing-to-other-ranks into a single buffer (to root)
  std::vector<long> keep_idx;
  keep_idx.reserve(static_cast<std::size_t>(_n_local));

  std::vector<double> out_buf;
  out_buf.reserve(static_cast<std::size_t>(_n_local) * STRIDE / 4);

  unsigned int out_count = 0;

  for (long i = 0; i < _n_local; ++i)
  {
    const long dr = da[i];

    if (dr < 0)
    {
      if (_drop_out_of_mesh)
        continue;
      keep_idx.push_back(i);
      continue;
    }

    if (static_cast<unsigned int>(dr) == static_cast<unsigned int>(_rank))
    {
      keep_idx.push_back(i);
      continue;
    }

    if (dr >= static_cast<long>(_n_procs))
    {
      if (_drop_out_of_mesh)
        continue;
      keep_idx.push_back(i);
      continue;
    }

    // send to root for routing
    out_buf.push_back(static_cast<double>(dr));
    out_buf.push_back(xa[i][0]); out_buf.push_back(xa[i][1]); out_buf.push_back(xa[i][2]);
    out_buf.push_back(va[i][0]); out_buf.push_back(va[i][1]); out_buf.push_back(va[i][2]);
    out_buf.push_back(oa[i][0]); out_buf.push_back(oa[i][1]); out_buf.push_back(oa[i][2]);
    out_buf.push_back(fa[i][0]); out_buf.push_back(fa[i][1]); out_buf.push_back(fa[i][2]);
    out_buf.push_back(ma[i][0]);

    ++out_count;
  }

  // Receive buffer from root (particles routed to this rank)
  unsigned int in_count = 0;
  std::vector<double> in_buf;

  const auto root = static_cast<libMesh::processor_id_type>(_root_rank);

  if (_rank != root)
  {
    // Always send count
    comm.send(root, out_count, tag_out_count);
    if (out_count > 0)
      comm.send(root, out_buf, tag_out_data);

    // Always receive count
    comm.receive(root, in_count, tag_in_count);
    if (in_count > 0)
    {
      in_buf.resize(static_cast<std::size_t>(in_count) * STRIDE);
      comm.receive(root, in_buf, tag_in_data);
    }
  }
  else
  {
    // Root receives outgoing from all ranks, routes by dest, then sends back to each rank.

    std::vector<std::vector<double>> route_buf(_n_procs);
    std::vector<unsigned int> route_count(_n_procs, 0);

    // Include root's own outgoing too (already packed)
    if (out_count > 0)
    {
      // unpack and push to destination buckets
      for (unsigned int k = 0; k < out_count; ++k)
      {
        const std::size_t base = static_cast<std::size_t>(k) * STRIDE;
        const int dest = static_cast<int>(out_buf[base + 0] + 0.5);
        if (dest < 0 || dest >= static_cast<int>(_n_procs))
          continue;
        auto & rb = route_buf[static_cast<std::size_t>(dest)];
        rb.insert(rb.end(), out_buf.begin() + base, out_buf.begin() + base + STRIDE);
        route_count[static_cast<std::size_t>(dest)] += 1u;
      }
    }

    // Receive from other ranks
    for (unsigned int src = 0; src < _n_procs; ++src)
    {
      if (src == _rank)
        continue;

      unsigned int csrc = 0;
      comm.receive(src, csrc, tag_out_count);

      std::vector<double> buf;
      if (csrc > 0)
      {
        buf.resize(static_cast<std::size_t>(csrc) * STRIDE);
        comm.receive(src, buf, tag_out_data);

        for (unsigned int k = 0; k < csrc; ++k)
        {
          const std::size_t base = static_cast<std::size_t>(k) * STRIDE;
          const int dest = static_cast<int>(buf[base + 0] + 0.5);
          if (dest < 0 || dest >= static_cast<int>(_n_procs))
            continue;
          auto & rb = route_buf[static_cast<std::size_t>(dest)];
          rb.insert(rb.end(), buf.begin() + base, buf.begin() + base + STRIDE);
          route_count[static_cast<std::size_t>(dest)] += 1u;
        }
      }
    }

    // Now send to each rank its routed particles
    for (unsigned int dest = 0; dest < _n_procs; ++dest)
    {
      if (dest == _rank)
      {
        in_count = route_count[dest];
        in_buf.swap(route_buf[dest]);
        continue;
      }

      comm.send(dest, route_count[dest], tag_in_count);
      if (route_count[dest] > 0)
        comm.send(dest, route_buf[dest], tag_in_data);
    }
  }

  // Build new local arrays = kept + received (strip dest field)
  const long n_recv = static_cast<long>(in_count);
  const long Nnew = static_cast<long>(keep_idx.size()) + n_recv;

  std::vector<double> x_new(static_cast<std::size_t>(Nnew) * 3, 0.0);
  std::vector<double> v_new(static_cast<std::size_t>(Nnew) * 3, 0.0);
  std::vector<double> o_new(static_cast<std::size_t>(Nnew) * 3, 0.0);
  std::vector<double> f_new(static_cast<std::size_t>(Nnew) * 3, 0.0);
  std::vector<double> m_new(static_cast<std::size_t>(Nnew) * 1, 0.0);

  long out = 0;

  for (const long i : keep_idx)
  {
    x_new[3*out+0]=xa[i][0]; x_new[3*out+1]=xa[i][1]; x_new[3*out+2]=xa[i][2];
    v_new[3*out+0]=va[i][0]; v_new[3*out+1]=va[i][1]; v_new[3*out+2]=va[i][2];
    o_new[3*out+0]=oa[i][0]; o_new[3*out+1]=oa[i][1]; o_new[3*out+2]=oa[i][2];
    f_new[3*out+0]=fa[i][0]; f_new[3*out+1]=fa[i][1]; f_new[3*out+2]=fa[i][2];
    m_new[out]=ma[i][0];
    ++out;
  }

  for (unsigned int k = 0; k < in_count; ++k)
  {
    const std::size_t base = static_cast<std::size_t>(k) * STRIDE;

    // base+0 is dest rank; ignore
    x_new[3*out+0]=in_buf[base+1]; x_new[3*out+1]=in_buf[base+2]; x_new[3*out+2]=in_buf[base+3];
    v_new[3*out+0]=in_buf[base+4]; v_new[3*out+1]=in_buf[base+5]; v_new[3*out+2]=in_buf[base+6];
    o_new[3*out+0]=in_buf[base+7]; o_new[3*out+1]=in_buf[base+8]; o_new[3*out+2]=in_buf[base+9];
    f_new[3*out+0]=in_buf[base+10]; f_new[3*out+1]=in_buf[base+11]; f_new[3*out+2]=in_buf[base+12];
    m_new[out]=in_buf[base+13];
    ++out;
  }

  if (out != Nnew)
    mooseError("TorchParticleCloudUserObject: migrateParticlesViaRoot internal size mismatch.");

  _n_local = Nnew;

  auto cpu = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);

  _x     = torch::from_blob(x_new.data(), {_n_local, 3}, cpu).clone().to(_device);
  _v_p   = torch::from_blob(v_new.data(), {_n_local, 3}, cpu).clone().to(_device);
  _omega = torch::from_blob(o_new.data(), {_n_local, 3}, cpu).clone().to(_device);
  _f_drag= torch::from_blob(f_new.data(), {_n_local, 3}, cpu).clone().to(_device);
  _mass  = torch::from_blob(m_new.data(), {_n_local, 1}, cpu).clone().to(_device);

  _cell_idx  = torch::full({_n_local}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));
  _dest_rank = torch::full({_n_local}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));
}

void
TorchParticleCloudUserObject::computeCellConcentrations()
{
  if (_local_elems.empty() || _n_local <= 0)
    return;

  const long Nc = static_cast<long>(_local_elems.size());
  if (!_cell_volume_t.defined() || _cell_volume_t.numel() != Nc)
    return;

  auto fopts = torch::TensorOptions().dtype(torch::kFloat64).device(_device);

  auto valid = _cell_idx.ge(0);
  const auto n_valid = valid.sum().item<int64_t>();
  if (n_valid <= 0)
    return;

  auto idx_valid  = _cell_idx.masked_select(valid).to(torch::kInt64);
  auto mass_valid = _mass.view({-1}).masked_select(valid) * _parcel_weight;
  auto num_valid  = torch::full({idx_valid.size(0)}, _parcel_weight, fopts);

  auto cell_mass = torch::zeros({Nc}, fopts);
  auto cell_num  = torch::zeros({Nc}, fopts);

  cell_mass.index_add_(0, idx_valid, mass_valid);
  cell_num.index_add_(0, idx_valid, num_valid);

  auto cell_mconc = cell_mass / (_cell_volume_t + 1e-30);
  auto cell_nconc = cell_num  / (_cell_volume_t + 1e-30);

  auto m_cpu = cell_mconc.to(torch::kCPU).contiguous();
  auto n_cpu = cell_nconc.to(torch::kCPU).contiguous();

  auto macc = m_cpu.accessor<double, 1>();
  auto nacc = n_cpu.accessor<double, 1>();

  for (long c = 0; c < Nc; ++c)
  {
    _cell_mass_conc[static_cast<std::size_t>(c)]   = static_cast<Real>(macc[c]);
    _cell_number_conc[static_cast<std::size_t>(c)] = static_cast<Real>(nacc[c]);
  }
}

void
TorchParticleCloudUserObject::chooseSubstepping(const Real dt_step, unsigned int & n_sub, Real & dt_sub) const
{
  n_sub = 1;
  dt_sub = dt_step;

  if (_max_cfl <= 0.0 || _max_substeps <= 1 || _min_cell_h <= 0.0)
    return;

  Real vmax = 0.0;
  if (_v_p.defined() && _v_p.numel() > 0)
    vmax = _v_p.norm(2, 1).max().item<double>();

  const Real dt_target = _max_cfl * _min_cell_h / (vmax + 1e-12);
  if (dt_target >= dt_step)
    return;

  const Real ratio = dt_step / std::max(dt_target, 1e-30);
  n_sub = std::min(_max_substeps, static_cast<unsigned int>(std::ceil(ratio)));
  n_sub = std::max(1u, n_sub);
  dt_sub = dt_step / static_cast<Real>(n_sub);
}

long
TorchParticleCloudUserObject::injectedTotalParcels(const Real time) const
{
  if (!_enable_injection || _injection_rate_parcels <= 0.0)
    return 0;

  const Real t0 = _injection_start_time;

  if (!std::isfinite(time) || time <= t0)
    return 0;

  Real t1 = time;
  if (_injection_end_time >= 0.0)
    t1 = std::min(t1, _injection_end_time);

  if (t1 <= t0)
    return 0;

  const Real s = t1 - t0;
  const Real n_real = _injection_rate_parcels * s;

  const long n = static_cast<long>(std::floor(n_real + 1e-12));
  return std::max(0L, n);
}

void
TorchParticleCloudUserObject::appendParticles(const long n_add,
                                             const Point & box_min,
                                             const Point & box_max,
                                             const RealVectorValue & vel0)
{
  if (n_add <= 0)
    return;

  // Ensure tensors exist (even if empty)
  auto opts = torch::TensorOptions().dtype(torch::kFloat64).device(_device);
  if (!_x.defined())
    _x = torch::zeros({0, 3}, opts);
  if (!_v_p.defined())
    _v_p = torch::zeros({0, 3}, opts);
  if (!_omega.defined())
    _omega = torch::zeros({0, 3}, opts);
  if (!_mass.defined())
    _mass = torch::zeros({0, 1}, opts);
  if (!_f_drag.defined())
    _f_drag = torch::zeros({0, 3}, opts);
  if (!_cell_idx.defined())
    _cell_idx = torch::full({0}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));
  if (!_dest_rank.defined())
    _dest_rank = torch::full({0}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));

  auto almost_equal = [](Real a, Real b) { return std::abs(a - b) < 1e-14; };

  // Deterministic per-rank, per-step host RNG (does not touch torch RNG state)
  const uint64_t seed =
      static_cast<uint64_t>(_rng_seed) ^
      (0x9e3779b97f4a7c15ULL + static_cast<uint64_t>(_rank) * 0xBF58476D1CE4E5B9ULL +
       static_cast<uint64_t>(_step_counter) * 0x94D049BB133111EBULL);
  std::mt19937_64 rng(seed);

  std::uniform_real_distribution<double> dist01(0.0, 1.0);

  std::vector<double> x_add(static_cast<std::size_t>(n_add) * 3, 0.0);
  std::vector<double> v_add(static_cast<std::size_t>(n_add) * 3, 0.0);
  std::vector<double> o_add(static_cast<std::size_t>(n_add) * 3, 0.0);
  std::vector<double> f_add(static_cast<std::size_t>(n_add) * 3, 0.0);
  std::vector<double> m_add(static_cast<std::size_t>(n_add) * 1, 0.0);

  const Real volume = (libMesh::pi / 6.0) * std::pow(_particle_diameter, 3);
  const Real m = _particle_density * volume;

  for (long i = 0; i < n_add; ++i)
  {
    for (int d = 0; d < 3; ++d)
    {
      const double lo = box_min(d);
      const double hi = box_max(d);

      // If box is degenerate in this direction, keep at lo
      const double r = dist01(rng);
      x_add[3 * i + d] = (almost_equal(lo, hi) ? lo : lo + (hi - lo) * r);
    }

    v_add[3 * i + 0] = vel0(0);
    v_add[3 * i + 1] = vel0(1);
    v_add[3 * i + 2] = vel0(2);

    // omega, f_drag already zero
    m_add[i] = m;
  }

  auto cpu = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
  auto x_t = torch::from_blob(x_add.data(), {n_add, 3}, cpu).clone().to(_device);
  auto v_t = torch::from_blob(v_add.data(), {n_add, 3}, cpu).clone().to(_device);
  auto o_t = torch::from_blob(o_add.data(), {n_add, 3}, cpu).clone().to(_device);
  auto f_t = torch::from_blob(f_add.data(), {n_add, 3}, cpu).clone().to(_device);
  auto m_t = torch::from_blob(m_add.data(), {n_add, 1}, cpu).clone().to(_device);

  auto idx_t = torch::full({n_add}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));
  auto dst_t = torch::full({n_add}, -1, torch::TensorOptions().dtype(torch::kInt64).device(_device));

  // Append
  _x     = torch::cat({_x, x_t}, 0);
  _v_p   = torch::cat({_v_p, v_t}, 0);
  _omega = torch::cat({_omega, o_t}, 0);
  _f_drag= torch::cat({_f_drag, f_t}, 0);
  _mass  = torch::cat({_mass, m_t}, 0);

  _cell_idx  = torch::cat({_cell_idx, idx_t}, 0);
  _dest_rank = torch::cat({_dest_rank, dst_t}, 0);

  _n_local += n_add;
}

void
TorchParticleCloudUserObject::injectParticles()
{
  if (!_enable_injection || _injection_rate_parcels <= 0.0)
    return;

  const Real t_prev = _t - _dt;

  const long n_prev = injectedTotalParcels(t_prev);
  const long n_now  = injectedTotalParcels(_t);
  const long n_add_global = std::max(0L, n_now - n_prev);

  if (n_add_global <= 0)
    return;

  const long base = n_add_global / static_cast<long>(_n_procs);
  const long rem  = n_add_global % static_cast<long>(_n_procs);
  const long n_add_local = base + ((static_cast<long>(_rank) < rem) ? 1 : 0);

  if (n_add_local <= 0)
    return;

  // Use injection box if provided; otherwise fall back to init_box (common use case)
  Point bmin = _injection_box_min;
  Point bmax = _injection_box_max;

  auto box_degenerate = [&](const Point & a, const Point & b)
  {
    for (int d = 0; d < 3; ++d)
      if (std::abs(a(d) - b(d)) > 1e-14)
        return false;
    return true;
  };

  if (box_degenerate(bmin, bmax) && !box_degenerate(_init_box_min, _init_box_max))
  {
    bmin = _init_box_min;
    bmax = _init_box_max;
  }

  if (box_degenerate(bmin, bmax))
    mooseError("TorchParticleCloudUserObject: injection is enabled but injection_box is degenerate (set injection_box_min/max or init_box_min/max).");

  RealVectorValue v0 = _injection_velocity;
  if (std::abs(v0(0)) < 1e-14 && std::abs(v0(1)) < 1e-14 && std::abs(v0(2)) < 1e-14 &&
      (std::abs(_init_velocity(0)) > 1e-14 || std::abs(_init_velocity(1)) > 1e-14 || std::abs(_init_velocity(2)) > 1e-14))
    v0 = _init_velocity;

  appendParticles(n_add_local, bmin, bmax, v0);
}

void
TorchParticleCloudUserObject::compressParticles(const torch::Tensor & keep_mask)
{
  if (_n_local <= 0)
    return;

  if (!keep_mask.defined() || keep_mask.numel() != _n_local)
    mooseError("TorchParticleCloudUserObject: compressParticles keep_mask size mismatch.");

  auto idx = torch::nonzero(keep_mask).view({-1});
  const long n_keep = idx.numel();

  if (n_keep == _n_local)
    return;

  auto ix = idx.to(_device);

  _x     = _x.index_select(0, ix);
  _v_p   = _v_p.index_select(0, ix);
  _omega = _omega.index_select(0, ix);
  _mass  = _mass.index_select(0, ix);
  _f_drag= _f_drag.index_select(0, ix);

  _cell_idx  = _cell_idx.index_select(0, ix);
  _dest_rank = _dest_rank.index_select(0, ix);

  _n_local = n_keep;
}

void
TorchParticleCloudUserObject::applyWallDeposition()
{
  _n_deposited_step = 0;

  if (_wall_deposition_velocity <= 0.0)
    return;

  const Real dt_step = _dt;
  if (!std::isfinite(dt_step) || dt_step <= 0.0)
    return;

  if (_n_local <= 0)
    return;

  const long n_cells = static_cast<long>(_local_elems.size());
  if (n_cells <= 0)
    return;

  if (!_cell_aw_over_v_t.defined() || _cell_aw_over_v_t.numel() != n_cells)
    return;

  // Only boundary-adjacent cells can deposit; if there are none, exit early
  if (_cell_aw_over_v_t.max().item<double>() <= 0.0)
    return;

  auto opts = torch::TensorOptions().dtype(torch::kFloat64).device(_device);

  // Per-cell effective deposition velocity
  torch::Tensor vd_cell = torch::full({n_cells}, _wall_deposition_velocity, opts);

  if (_include_settling_in_vd)
  {
    const Real gmag = std::sqrt(_gravity(0) * _gravity(0) + _gravity(1) * _gravity(1) + _gravity(2) * _gravity(2));
    const Real d = _particle_diameter;

    std::vector<double> vs_host(static_cast<std::size_t>(n_cells), 0.0);
    for (long c = 0; c < n_cells; ++c)
    {
      const Real rho_f = std::max(_cell_rho[static_cast<std::size_t>(c)], 1e-30);
      const Real mu_f  = std::max(_cell_mu[static_cast<std::size_t>(c)], _mu_min);

      // Stokes settling speed (no slip correction); clamp to [0, v_max]
      Real vs = ((_particle_density - rho_f) * gmag * d * d) / (18.0 * mu_f + 1e-30);
      vs = std::max(0.0, std::min(vs, _v_max));
      vs_host[static_cast<std::size_t>(c)] = static_cast<double>(vs);
    }

    auto cpu = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
    auto vs_t = torch::from_blob(vs_host.data(), {n_cells}, cpu).clone().to(_device);
    vd_cell = vd_cell + vs_t;
  }

  // k_dep = Vd * (A_wall/V_cell) [1/s]
  auto k_cell = vd_cell * _cell_aw_over_v_t;
  auto valid = _cell_idx.ge(0);

  auto idx_safe = _cell_idx.clamp(0, n_cells - 1);
  auto k_p = k_cell.index_select(0, idx_safe).view({-1, 1});

  auto p = 1.0 - torch::exp(-k_p * dt_step);
  p = torch::clamp(p, 0.0, 1.0) * valid.to(torch::kFloat64).view({-1, 1});

  if (p.max().item<double>() <= 0.0)
    return;

  // Stochastic deposition
  auto r = torch::rand({_n_local, 1}, opts);
  auto dep = r.lt(p).view({-1});

  const long long n_dep = dep.sum().item<int64_t>();
  if (n_dep <= 0)
    return;

  _n_deposited_step = static_cast<unsigned long long>(n_dep);
  _n_deposited_total += _n_deposited_step;

  auto keep = dep.logical_not();
  compressParticles(keep);
}


void
TorchParticleCloudUserObject::advanceParticlesExplicit()
{
  const Real dt_step = _dt;

  if (!std::isfinite(dt_step) || dt_step <= 0.0)
  {
    if (_f_drag.defined())
      _f_drag.zero_();
    return;
  }

  if (_n_local <= 0)
    return;

  const long n_cells = static_cast<long>(_local_elems.size());
  if (n_cells <= 0)
  {
    _f_drag.zero_();
    return;
  }

  unsigned int n_sub = 1;
  Real dt_sub = dt_step;
  chooseSubstepping(dt_step, n_sub, dt_sub);

  auto opts = torch::TensorOptions().dtype(torch::kFloat64).device(_device);

  std::vector<double> vel_host(3 * _local_elems.size(), 0.0);
  std::vector<double> rho_host(_local_elems.size(), 0.0);
  std::vector<double> mu_host(_local_elems.size(), 0.0);

  for (unsigned int i = 0; i < _local_elems.size(); ++i)
  {
    vel_host[3*i + 0] = _cell_vel[i](0);
    vel_host[3*i + 1] = _cell_vel[i](1);
    vel_host[3*i + 2] = _cell_vel[i](2);

    rho_host[i] = std::max(_cell_rho[i], 1e-30);
    mu_host[i]  = std::max(_cell_mu[i],  _mu_min);
  }

  auto cell_vel = torch::from_blob(vel_host.data(), {n_cells, 3}, torch::kFloat64).clone().to(_device);
  auto cell_rho = torch::from_blob(rho_host.data(), {n_cells, 1}, torch::kFloat64).clone().to(_device);
  auto cell_mu  = torch::from_blob(mu_host.data(),  {n_cells, 1}, torch::kFloat64).clone().to(_device);

  auto valid = _cell_idx.ge(0);
  auto idx_safe = _cell_idx.clamp(0, n_cells - 1);
  auto valid_f = valid.to(torch::kFloat64).view({-1, 1});

  auto u_f   = cell_vel.index_select(0, idx_safe);
  auto rho_f = cell_rho.index_select(0, idx_safe);
  auto mu_f  = cell_mu.index_select(0, idx_safe);

  const double d = _particle_diameter;
  const double A = libMesh::pi * 0.25 * d * d;
  auto A_t = torch::full({_n_local, 1}, A, opts);
  auto d_t = torch::full({_n_local, 1}, d, opts);

  auto g = torch::tensor({_gravity(0), _gravity(1), _gravity(2)}, opts).view({1, 3});

  auto f_drag_acc = torch::zeros_like(_f_drag);

  const double kmin = 1e-12;

  for (unsigned int s = 0; s < n_sub; ++s)
  {
    auto rel = u_f - _v_p;
    auto rel_mag = rel.norm(2, 1).view({-1, 1});

    auto Re = (rho_f * rel_mag * d_t) / (mu_f + 1e-30);
    Re = torch::clamp(Re, 0.0, _re_max);

    auto Cd_stokes = 24.0 / (Re + 1e-30);
    auto Cd = Cd_stokes * (1.0 + 0.15 * torch::pow(Re, 0.687));
    Cd = torch::where(Re > 1000.0, torch::full_like(Cd, 0.44), Cd);
    Cd = torch::clamp(Cd, 0.0, _cd_max);

    auto K = (0.5 * Cd * rho_f * A_t * rel_mag) / (_mass + 1e-30);
    K = torch::clamp(K, 0.0, _k_max);

    auto a_body = (1.0 - (rho_f / _particle_density)) * g;

    auto Kdt = torch::clamp(K * dt_sub, 0.0, 700.0);
    auto ef = torch::exp(-Kdt);
    auto phi = torch::where(K > kmin, (1.0 - ef) / (K + 1e-30), torch::full_like(K, dt_sub));

    auto v_new = u_f + (_v_p - u_f) * ef + a_body * phi;

    auto dt_t = torch::full_like(K, dt_sub);
    auto term3 = torch::where(K > kmin, a_body * (dt_t - phi) / (K + 1e-30), 0.5 * a_body * dt_t * dt_t);
    auto x_new = _x + u_f * dt_sub + (_v_p - u_f) * phi + term3;

    _v_p = _v_p + (v_new - _v_p) * valid_f;
    _x   = _x   + (x_new - _x)   * valid_f;

    auto speed = _v_p.norm(2, 1).view({-1, 1});
    auto too_fast = speed.gt(_v_max);
    _v_p = torch::where(too_fast, _v_p * (_v_max / (speed + 1e-30)), _v_p);

    auto rel_new = u_f - _v_p;
    auto f_drag = (_mass * K) * rel_new;

    f_drag_acc = f_drag_acc + f_drag * valid_f;
  }

  _f_drag = f_drag_acc / static_cast<double>(n_sub);
  _f_drag = torch::where(valid.view({-1, 1}), _f_drag, torch::zeros_like(_f_drag));
}

#endif
