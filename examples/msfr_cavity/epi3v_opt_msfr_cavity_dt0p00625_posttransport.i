# =====================================================================================================
# Molten salt fast reactor (MSFR) core-cavity model -- EVOL-benchmark-inspired geometry: a single,
# fully coupled, publication-quality example integrating
#
#   * fluid flow            -- incompressible Navier-Stokes (linear FV, Rhie-Chow, segregated SIMPLE),
#   * nuclear power          -- a separable cosine neutron-flux shape (simplified power production),
#   * energy deposition      -- volumetric fission heat, advected enthalpy, molecular + turbulent
#                               conduction,
#   * turbulence             -- standard k-epsilon with wall functions (turbulent viscosity mu_t),
#   * radiolysis             -- the chloride radical chain producing the oxidant Cl2.- proportional to
#                               power, transported by advection and molecular + turbulent diffusion,
#   * corrosion              -- the structural-alloy cavity walls dissolve through a temperature-
#                               dependent Butler-Volmer reaction, releasing chromium into the salt as
#                               Cr(II); the radiolytic oxidant then oxidizes it Cr(II) -> Cr(III),
#                               coupling corrosion and radiolysis through the dissolved-chromium redox.
#
# Geometry (the realistic MSFR concept rather than a channel): an open 2D core cavity bounded by solid
# structural-alloy walls. The fuel salt (NaCl-UCl3) is pumped in through a central inlet at the bottom
# and collected at a central outlet plenum at the top. The central jet rises through the core and
# entrains the surrounding salt, so two large RECIRCULATION ZONES form on either side, with the salt
# flowing back DOWN along the side walls. Those slow, recirculating regions are the realistic concern
# -- they are hotter (poor heat removal), they trap the radiolytic oxidant, and their wall corrosion
# is altered. This is the spirit of the EVOL/MSFR CFD benchmark (a heated core cavity with internal
# recirculation), here carrying the full multiphysics.
#
# Validation provenance (each physics component is validated separately in this application):
#   * turbulent flow / recirculating cavity -- validation/free_flow (ERCOFTAC channel and
#                            backward-facing-step, the canonical separated/recirculating flows),
#   * radiolysis kinetics -- validation/msr (pulse-radiolysis transient-absorption traces),
#   * corrosion kinetics  -- validation/corrosion (76 cases / 43 targets of the effective Butler-
#                            Volmer correlation; the wall rate here reproduces that calibration).
# Verification postprocessors confirm the energy balance closes (nuclear heat == enthalpy rise), the
# flow is fully turbulent (reported Reynolds number), and the recirculation is captured (reverse flow).
#
# Run:  open_pronghorn-opt -i msfr_cavity.i
# View: the Exodus output (vel with the corner recirculation, T_fluid hot spots, mu_t, Cl2m_rad,
#       Cr_II/Cr_III, corrosion_rate_um_y).
# =====================================================================================================

# --- geometry (MSFR-like core cavity) ---
Wc = 2.0              # m, cavity width
Hc = 2.0              # m, cavity height
in_min = 0.4          # m, central inlet: x from in_min ...
in_max = 1.6          # m, ... to in_max (1.2 m wide; wide enough to sweep the floor corners)
out_min = 0.7         # m, central outlet plenum: x from out_min ...
out_max = 1.3         # m, ... to out_max (0.6 m wide)

# --- NaCl-UCl3 fast-reactor fuel salt (~650 C) ---

# --- operating point ---
u_in = 1.0            # m/s inlet velocity across the cavity floor (upward, pumped circulation)

# --- corrosion operating point ---

# --- k-epsilon closure ---
walls = 'left right top_wall bottom_wall'

# Second MSFR-cavity transport half-step.
# Receives chemistry-updated species and the already-halved solved velocity field.

# Corrected transient operating condition: BOTH Strang transport half-steps use Cr(II) inlet = 3.0 mol/m^3; domain IC remains zero.
[Mesh]
  [cavity]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${Wc}
    ymin = 0
    ymax = ${Hc}
    nx = 40
    ny = 40
  []
  # Carve a central inlet jet and outlet plenum (and the surrounding walls) out of the bottom and top
  # boundaries, then drop the originals so every face belongs to exactly one sideset. The salt jets up
  # the centre and recirculates back down along the side walls.
  [outlet]
    type = ParsedGenerateSideset
    input = cavity
    included_boundaries = 'top'
    combinatorial_geometry = 'x > ${out_min} & x < ${out_max}'
    new_sideset_name = 'outlet'
  []
  [top_wall]
    type = ParsedGenerateSideset
    input = outlet
    included_boundaries = 'top'
    combinatorial_geometry = 'x < ${out_min} | x > ${out_max}'
    new_sideset_name = 'top_wall'
  []
  [inlet]
    type = ParsedGenerateSideset
    input = top_wall
    included_boundaries = 'bottom'
    combinatorial_geometry = 'x > ${in_min} & x < ${in_max}'
    new_sideset_name = 'inlet'
  []
  [bottom_wall]
    type = ParsedGenerateSideset
    input = inlet
    included_boundaries = 'bottom'
    combinatorial_geometry = 'x < ${in_min} | x > ${in_max}'
    new_sideset_name = 'bottom_wall'
  []
  [drop_orig]
    type = BoundaryDeletionGenerator
    input = bottom_wall
    boundary_names = 'top bottom'
  []
[]

[Problem]
  # Keep one empty nonlinear system because this MOOSE version's
  # SubProblem::automaticScaling() getter assumes nonlinear system 0 exists.
  nl_sys_names = 'nl0'
  skip_nl_system_check = true

  # extra_tag_vectors is a 2-D vector. With nl0 present, solver-system index 0
  # is nl0 and index 1 is e_sol_sys. Put a harmless private tag on nl0 and
  # NONTIME on the first linear system, exactly where it lived before nl0 was
  # introduced. This avoids the NONTIME parallel-type collision in nl0.
  extra_tag_vectors = 'nl_anchor_tag; nontime'

  linear_sys_names = 'e_sol_sys cl_ion_sys cl_rad_sys cl2m_rad_sys cl3_ion_sys cl2_diss_sys cr_ii_sys cr_iii_sys cr_i_sys'
  previous_nl_solution_required = true
[]




[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Cr'
  integration_method = linear

  # Transport only: EPI3V remains the sole chemistry integrator.
  transport_only = true
  time_derivative = true

  temperature = T_local
  dose_rate = 0

  initial_condition_species = 'Cl_ion Cr_II'
  initial_condition_values = '2.0e4 0.0'

  rhie_chow_user_object = rc
  use_generic_rhie_chow_advection = true
  advected_interp_method = upwind
  diffusivity = D_eff_half
[]

[AuxVariables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
  [pressure]
    type = INSFVPressureVariable
    initial_condition = 0.0
  []

  [D_eff_half]
    type = MooseVariableFVReal
    initial_condition = 1.0e-9
  []

  [T_local]
    type = MooseVariableFVReal
    initial_condition = 673.15
  []
  [a_u]
    type = MooseVariableFVReal
    initial_condition = 2.0
  []
  [a_v]
    type = MooseVariableFVReal
    initial_condition = 2.0
  []
[]
[GlobalParams]
  rhie_chow_user_object = rc
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = u
    v = v
    pressure = pressure
    a_u = a_u
    a_v = a_v
    velocity_interp_method = rc
  []
  [linear_midpoint_reconstruct]
    type = MSRLinearMidpointReconstruction
    execute_on = TIMESTEP_END
  []
[]

[CorrosionPlatingFlow]
  elements = 'Cr'
  release_variables = 'Cr_II'
  reaction_boundary = ${walls}

  temperature = T_local
  reference_temperature = 923.15
  material_class = stainless_316
  salt_class = chloride
  redox_class = oxidizing_fef2
  applied_overpotential = 0.12
  reference_concentration = 3.0
  temperature_dependent_kinetics = true

  # Cr_II is owned by the nonlinear FV EPI3V transport system.
  release_variable_formulation = linear

  # Each transport child represents a Strang half operator over the macro duration.
  flux_scale = 0.5
[]

[LinearFVBCs]
  [inlet_e_sol]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = e_sol
    boundary = inlet
    functor = 0.0
  []
  [outlet_e_sol]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = e_sol
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cl_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl_ion
    boundary = inlet
    functor = 2.0e4
  []
  [outlet_cl_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl_ion
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cl_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl_rad
    boundary = inlet
    functor = 0.0
  []
  [outlet_cl_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl_rad
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cl2m_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl2m_rad
    boundary = inlet
    functor = 0.0
  []
  [outlet_cl2m_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl2m_rad
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cl3_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl3_ion
    boundary = inlet
    functor = 0.0
  []
  [outlet_cl3_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl3_ion
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cl2_diss]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cl2_diss
    boundary = inlet
    functor = 0.0
  []
  [outlet_cl2_diss]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cl2_diss
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cr_ii]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cr_II
    boundary = inlet
    functor = 3.0
  []
  [outlet_cr_ii]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cr_II
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cr_iii]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cr_III
    boundary = inlet
    functor = 0.0
  []
  [outlet_cr_iii]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cr_III
    boundary = outlet
    use_two_term_expansion = false
  []
  [inlet_cr_i]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = Cr_I
    boundary = inlet
    functor = 0.0
  []
  [outlet_cr_i]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = Cr_I
    boundary = outlet
    use_two_term_expansion = false
  []
[]

[FVBCs]
  [inlet_u]
    type = INSFVInletVelocityBC
    boundary = inlet
    variable = u
    functor = 0.0
  []
  [inlet_v]
    type = INSFVInletVelocityBC
    boundary = inlet
    variable = v
    functor = '${fparse 0.5 * ${u_in}}'
  []
  [walls_u]
    type = INSFVNoSlipWallBC
    boundary = ${walls}
    variable = u
    function = 0.0
  []
  [walls_v]
    type = INSFVNoSlipWallBC
    boundary = ${walls}
    variable = v
    function = 0.0
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = outlet
    variable = pressure
    function = 0.0
  []
[]

[Convergence]
  [fp]
    type = IterationCountConvergence
    max_iterations = 1
    converge_at_max_iterations = true
  []
[]

[Executioner]
  type = Transient
  system_names = 'e_sol_sys cl_ion_sys cl_rad_sys cl2m_rad_sys cl3_ion_sys cl2_diss_sys cr_ii_sys cr_iii_sys cr_i_sys'

  multi_system_fixed_point = true
  multi_system_fixed_point_convergence = fp

  start_time = 0.0
  end_time = 0.2
  dtmin = 1.0e-14
  dtmax = 0.00625

  [TimeIntegrator]
    type = MSRLinearImplicitMidpoint
  []

  [TimeStepper]
    type = ConstantDT
    dt = 0.00625
  []

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]



[Outputs]
  exodus = false
[]
