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

# -------------------------------------------------------------------------------------
# EPI3V transport app on the exact MSFR cavity mesh.
# Receives the solved physical velocity field as u_full/v_full.
# u/v are one-half of that field so advancing over h represents T(h/2).
# Chemistry receives actual MSFR T_fluid and the original spatial MSFR dose shape.
# -------------------------------------------------------------------------------------

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

[Functions]

  # Exact spatial radiolytic-dose definition from the original msfr_cavity.i.
  [epi3v_flux_shape]
    type = ParsedFunction
    expression = 'cos(pi*(x-0.5*${Wc})/(1.4*${Wc})) * cos(pi*(y-0.5*${Hc})/(1.4*${Hc}))'
  []
  [epi3v_dose_field]
    type = ParsedFunction
    expression = '5.0e5 * f'
    symbol_names = 'f'
    symbol_values = 'epi3v_flux_shape'
  []
[]

[AuxVariables]
  [u_full]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
  [v_full]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
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
  [T_local]
    type = MooseVariableFVReal
    initial_condition = 673.15
  []

  [dose_local]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []

  [mu_t_full]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
  [D_eff_half]
    type = MooseVariableFVReal
    initial_condition = 1.0e-9
  []
  [a_u_full]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
  [a_v_full]
    type = MooseVariableFVReal
    initial_condition = 1.0
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
[AuxKernels]
  [half_u]
    type = ParsedAux
    variable = u
    coupled_variables = 'u_full'
    expression = '0.5*u_full'
    execute_on = TIMESTEP_BEGIN
  []
  [half_v]
    type = ParsedAux
    variable = v
    coupled_variables = 'v_full'
    expression = '0.5*v_full'
    execute_on = TIMESTEP_BEGIN
  []

  [dose_local_from_msfr_shape]
    type = FunctionAux
    variable = dose_local
    function = epi3v_dose_field
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []

  # MSFR effective species diffusivity:
  # D_eff = D_mol + mu_t/(rho*Sc_t)
  # The factor 0.5 represents the Strang transport half-step.
  [compute_D_eff_half]
    type = ParsedAux
    variable = D_eff_half
    coupled_variables = 'mu_t_full'
    expression = '0.5 * (2.0e-9 + mu_t_full/(3300.0*0.7))'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  # For the Strang half operator:
  # u_half = 0.5*u_full and D_RC,half = V/(2*A), so a_half = 2*A.
  [scale_a_u_for_half_rc]
    type = ParsedAux
    variable = a_u
    coupled_variables = 'a_u_full'
    expression = '2.0*a_u_full'
    execute_on = TIMESTEP_BEGIN
  []
  [scale_a_v_for_half_rc]
    type = ParsedAux
    variable = a_v
    coupled_variables = 'a_v_full'
    expression = '2.0*a_v_full'
    execute_on = TIMESTEP_BEGIN
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

[MultiApps]
  [chemistry]
    type = TransientMultiApp
    sub_cycling = true
    input_files = epi3v_opt_msfr_cavity_dt0p00625_chemistry.i
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [species_to_chemistry]
    type = MultiAppCopyTransfer
    to_multi_app = chemistry
    source_variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
    variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
  []
  [velocity_to_chemistry]
    type = MultiAppCopyTransfer
    to_multi_app = chemistry
    source_variable = 'u v'
    variable = 'u v'
  []
  [temperature_to_chemistry]
    type = MultiAppCopyTransfer
    to_multi_app = chemistry
    source_variable = T_local
    variable = T_local
  []
  [dose_to_chemistry]
    type = MultiAppCopyTransfer
    to_multi_app = chemistry
    source_variable = dose_local
    variable = dose_local
  []
  [species_from_chemistry]
    type = MultiAppCopyTransfer
    from_multi_app = chemistry
    source_variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
    variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
  []

  [D_eff_half_to_chemistry]
    type = MultiAppCopyTransfer
    to_multi_app = chemistry
    source_variable = D_eff_half
    variable = D_eff_half
  []
  [rc_fields_to_chemistry]
    type = MultiAppCopyTransfer
    to_multi_app = chemistry
    source_variable = 'pressure a_u a_v'
    variable = 'pressure a_u a_v'
  []
[]

[VectorPostprocessors]
  [child_half_reconstructed_rc_flux]
    type = MSRRhieChowFaceFluxVPP
    rhie_chow_user_object = rc
    internal_only = true
    execute_on = MULTIAPP_FIXED_POINT_END
  []
[]

[Postprocessors]
  [u_full_min]
    type = ElementExtremeValue
    variable = u_full
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [u_full_max]
    type = ElementExtremeValue
    variable = u_full
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [v_full_min]
    type = ElementExtremeValue
    variable = v_full
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [v_full_max]
    type = ElementExtremeValue
    variable = v_full
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [u_half_min]
    type = ElementExtremeValue
    variable = u
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [u_half_max]
    type = ElementExtremeValue
    variable = u
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [v_half_min]
    type = ElementExtremeValue
    variable = v
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [v_half_max]
    type = ElementExtremeValue
    variable = v
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []

  [T_local_min]
    type = ElementExtremeValue
    variable = T_local
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [T_local_max]
    type = ElementExtremeValue
    variable = T_local
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []

  [Cr_II_avg]
    type = ElementAverageValue
    variable = Cr_II
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_II_min]
    type = ElementExtremeValue
    variable = Cr_II
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_II_max]
    type = ElementExtremeValue
    variable = Cr_II
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_III_avg]
    type = ElementAverageValue
    variable = Cr_III
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_III_min]
    type = ElementExtremeValue
    variable = Cr_III
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_III_max]
    type = ElementExtremeValue
    variable = Cr_III
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_I_avg]
    type = ElementAverageValue
    variable = Cr_I
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_I_min]
    type = ElementExtremeValue
    variable = Cr_I
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [Cr_I_max]
    type = ElementExtremeValue
    variable = Cr_I
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [total_Cr_avg]
    type = ParsedPostprocessor
    expression = 'Cr_II_avg + Cr_III_avg + Cr_I_avg'
    pp_names = 'Cr_II_avg Cr_III_avg Cr_I_avg'
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []

  [dose_local_min]
    type = ElementExtremeValue
    variable = dose_local
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [dose_local_max]
    type = ElementExtremeValue
    variable = dose_local
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []

  [mu_t_full_min]
    type = ElementExtremeValue
    variable = mu_t_full
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [mu_t_full_max]
    type = ElementExtremeValue
    variable = mu_t_full
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [D_eff_half_min]
    type = ElementExtremeValue
    variable = D_eff_half
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [D_eff_half_max]
    type = ElementExtremeValue
    variable = D_eff_half
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []

  [a_u_half_min]
    type = ElementExtremeValue
    variable = a_u
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [a_u_half_max]
    type = ElementExtremeValue
    variable = a_u
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [a_v_half_min]
    type = ElementExtremeValue
    variable = a_v
    value_type = min
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
  []
  [a_v_half_max]
    type = ElementExtremeValue
    variable = a_v
    value_type = max
    execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
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
  csv = true
  exodus = true
  execute_on = 'INITIAL MULTIAPP_FIXED_POINT_END'
[]
