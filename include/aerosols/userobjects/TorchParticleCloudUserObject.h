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

// MOOSE
#include "GeneralUserObject.h"
#include "ADFunctorInterface.h"
#include "MooseFunctor.h"
#include "MooseTypes.h"

// libMesh
#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/parallel.h"
#include "libmesh/point.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/vector_value.h"

// torch
#include <torch/torch.h>

// STL
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

/**
 * TorchParticleCloudUserObject (MPI-capable, debug-stable)
 *
 * This version adds an optional "root-serialized" migration mode intended as a robust
 * debugging baseline when all-to-all particle exchange deadlocks.
 *
 * - If migration_via_root=true (recommended for now), particles that leave a rank are sent to root_rank,
 *   and root forwards them to the destination rank. This serializes MPI traffic and avoids pairwise deadlocks.
 * - Optional MPI barriers can be enabled to force strict stage synchronization for debugging.
 */
class TorchParticleCloudUserObject : public GeneralUserObject, public ADFunctorInterface
{
public:
  static InputParameters validParams();
  TorchParticleCloudUserObject(const InputParameters & parameters);

  void initialize() override;
  void execute() override;
  void finalize() override;

  const libMesh::VectorValue<Real> & cellForce(const libMesh::Elem * elem) const;
  libMesh::VectorValue<Real> cellForceDensity(const libMesh::Elem * elem) const;

  Real cellNumberConcentration(const libMesh::Elem * elem) const;
  Real cellMassConcentration(const libMesh::Elem * elem) const;

protected:
  void buildCellCache();
  void sampleFluidToCells();

  void ensureProcessorBoundingBoxes();
  libMesh::processor_id_type bboxRoute(const libMesh::Point & p) const;

  void mapParticlesToCells();

  // Migration (MPI)
  void migrateParticles();
  void migrateParticlesViaRoot(); // serialized through root_rank

  void computeCellConcentrations();
  void chooseSubstepping(const Real dt_step, unsigned int & n_sub, Real & dt_sub) const;
  void advanceParticlesExplicit();


  // Injection (optional)
  long injectedTotalParcels(const Real time) const;
  void injectParticles();
  void appendParticles(const long n_add,
                     const Point & box_min,
                     const Point & box_max,
                     const RealVectorValue & vel0);


  // Wall deposition (optional)
  void applyWallDeposition();
  void compressParticles(const torch::Tensor & keep_mask);

  // -----------------------------
  // Parameters
  // -----------------------------
  const unsigned int _mesh_dim;

  const unsigned int _num_particles_initial;
  const Real _particle_density;
  const Real _particle_diameter;
  const Real _parcel_weight;

  const Point _init_box_min;
  const Point _init_box_max;
  const RealVectorValue _init_velocity;


  // Injection (optional)
  const bool _enable_injection;
  const Real _injection_rate_parcels; // GLOBAL computational parcels per second
  const Real _injection_start_time;
  const Real _injection_end_time; // <0 => no end time (inject forever)
  const Point _injection_box_min;
  const Point _injection_box_max;
  const RealVectorValue _injection_velocity;


  // Wall deposition (optional)
  const Real _wall_deposition_velocity; // [m/s] effective deposition velocity
  const bool _include_settling_in_vd;

  const RealVectorValue _gravity;

  const unsigned int _max_substeps;
  const Real _max_cfl;

  const unsigned int _rng_seed;
  const bool _debug;

  // MPI behavior
  const bool _drop_out_of_mesh;
  const bool _bbox_fallback;
  const int _migration_tag_base;

  // Debugging synchronization / serialized migration
  const bool _mpi_barriers;
  const bool _migration_via_root;
  const unsigned int _root_rank;

  // Drag safety clamps
  const Real _mu_min;
  const Real _re_max;
  const Real _cd_max;
  const Real _k_max;
  const Real _v_max;

  // -----------------------------
  // Fluid functors
  // -----------------------------
  const Moose::Functor<ADReal> & _u;
  const Moose::Functor<ADReal> * _v;
  const Moose::Functor<ADReal> * _w;
  const Moose::Functor<ADReal> & _rho;
  const Moose::Functor<ADReal> & _mu;

  // -----------------------------
  // MPI / rank info
  // -----------------------------
  libMesh::processor_id_type _rank = 0;
  libMesh::processor_id_type _n_procs = 1;

  // Bookkeeping (no MPI needed)
  unsigned long long _step_counter = 0;
  mutable unsigned long long _n_deposited_total = 0;
  mutable unsigned long long _n_deposited_step = 0;
  mutable unsigned long long _n_dropped_total = 0;
  mutable unsigned long long _n_dropped_step = 0;

  bool _proc_bbox_built = false;
  std::vector<libMesh::Point> _proc_bbox_min;
  std::vector<libMesh::Point> _proc_bbox_max;

  // -----------------------------
  // Cell cache (host, local elems)
  // -----------------------------
  bool _cache_built = false;
  Real _min_cell_h = std::numeric_limits<Real>::max();

  std::vector<const libMesh::Elem *> _local_elems;
  std::vector<Point> _cell_centroid;
  std::vector<Real> _cell_volume;
  std::vector<Real> _cell_wall_area; // sum of boundary side measure (area in 3D, length in 2D)

  std::vector<RealVectorValue> _cell_vel;
  std::vector<Real> _cell_rho;
  std::vector<Real> _cell_mu;

  std::unordered_map<dof_id_type, unsigned int> _elem_id_to_local;
  std::unique_ptr<libMesh::PointLocatorBase> _point_locator;

  std::vector<Real> _cell_number_conc;
  std::vector<Real> _cell_mass_conc;

  torch::Tensor _cell_volume_t; // [Nc] float64 on device
  torch::Tensor _cell_aw_over_v_t; // [Nc] (A_wall/V_cell) on device

  // -----------------------------
  // Particle state (torch, local only)
  // -----------------------------
  torch::Device _device;

  long _n_local = 0; // local particle count on this rank

  torch::Tensor _x;       // [Nlocal,3]
  torch::Tensor _v_p;     // [Nlocal,3]
  torch::Tensor _omega;   // [Nlocal,3]
  torch::Tensor _mass;    // [Nlocal,1]

  torch::Tensor _cell_idx;   // [Nlocal] local cell index or -1
  torch::Tensor _dest_rank;  // [Nlocal] owning rank or -1

  torch::Tensor _f_drag;  // [Nlocal,3]

  // Outputs
  std::unordered_map<dof_id_type, libMesh::VectorValue<Real>> _cell_forces;
  static const libMesh::VectorValue<Real> _zero_vec;

  mutable Real _last_debug_time = std::numeric_limits<Real>::quiet_NaN();
};

#endif
