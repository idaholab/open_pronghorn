#===============================================================================
# Chen test – air + aerosol Lagrangian particles (torch)
#===============================================================================

# --- Fluid properties: AIR (approx. at 293 K)  ### CHANGED
rho = 1.2        # kg/m^3
mu  = 1.8e-5     # Pa·s
g   = '0 0 -9.81'

### k-epsilon Closure Parameters ###
sigma_k   = 0.7194
sigma_eps = 0.7194
C_mu      = 0.0845
C1_eps    = 1.42
C2_eps    = 1.68
bulk_wall_treatment = false
wall_treatment      = 'neq'

# --- Aerosol properties (Chen / Gray case)
d_p      = 1e-5       # particle diameter [m] = 1 micrometer
rho_p    = 1400.0     # particle density [kg/m^3]

# --- Inlet turbulence settings (from Gray’s Chen case)
U_in = 0.225
Re_dh = ${fparse rho*0.04*U_in/mu}
I_inlet = ${fparse 0.16*Re_dh^(-1./8.)} #0.02   # turbulence intensity
L_in = ${fparse 0.04*0.07}
k_in  = ${fparse 1.5 * (U_in * I_inlet)^2}
eps_in = ${fparse C_mu^(3/4) * k_in^(3/2) / L_in}
# For now, we keep the inlet TKE/TKED BCs as simple constants (see BCs).

#===============================================================================
# MESH
#===============================================================================
[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim  = 3
    dx   = '0.8'
    dy   = '0.18 0.04 0.18'
    dz   = '0.02 0.04 0.28 0.04 0.02'
    ix   = '30'
    iy   = '8 4 8'
    iz   = '5 3 10 3 5'
    subdomain_id = '0 0 0
                    0 1 0
                    0 0 0
                    0 2 0
                    0 0 0
                    '
  []
  [rename_boundaries]
    type = RenameBoundaryGenerator
    input = cmg
    old_boundary = 'back front top  bottom left right'
    new_boundary = 'wall wall  wall wall   wall wall'
  []
  [create_inlet]
    type = ParsedGenerateSideset
    input = rename_boundaries
    combinatorial_geometry = 'x < 1e-5'
    new_sideset_name = 'inlet'
    fixed_normal = true
    normal = '-1 0 0'
    replace = true
    included_subdomains = '2'
  []
  [create_outlet]
    type = ParsedGenerateSideset
    input = create_inlet
    combinatorial_geometry = 'x > (0.4 - 1e-5)'
    new_sideset_name = 'outlet'
    fixed_normal = true
    normal = '1 0 0'
    replace = true
    included_subdomains = '1'
  []
  [rename_blocks]
    type = RenameBlockGenerator
    input = create_outlet
    old_block = '0 1 2'
    new_block = '0 0 0'
  []
[]

#===============================================================================
# GLOBAL PARAMS (needed by FV kernels)
#===============================================================================
[GlobalParams]
  rhie_chow_user_object = 'RhieChowMassFlux'
  advected_interp_method = 'upwind'

  pressure = pressure

  u = vel_x
  v = vel_y
  w = vel_z

  rho = ${rho}
[]

#===============================================================================
# PROBLEM
#===============================================================================
[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system
                      TKE_system TKED_system'
  previous_nl_solution_required = true
  error_on_jacobian_nonzero_reallocation = false
  kernel_coverage_check = false
[]

[UserObjects]
  [RhieChowMassFlux]
    type = RhieChowMassFlux
    pressure = pressure
    p_diffusion_kernel = p_diffusion
  []
  [particle_cloud]
    type = TorchParticleCloudUserObject
    execute_on = 'initial timestep_end'

    # Fluid functors (FV variables / functors)
    u = vel_x
    v = vel_y
    w = vel_z
    rho = rho
    mu  = mu

    # Particle properties
    # parcel_weight = 1.0
    num_particles = 0 # initial particles

    # Injection: either physical #/s...
    injection_rate = 1e9
    parcel_weight = 1e6   # => 1e3 parcels/s
    injection_start_time = 0.0
    injection_end_time   = 100.0

    particle_density  = ${rho_p}
    particle_diameter = ${d_p}

    # Seed particles in a thin slab just inside the inlet
    # (chosen to match the inlet block around y=[0.18,0.22], z=[0.06,0.34])
    init_box_min = '0.00 0.18 0.34'
    init_box_max = '0.02 0.22 0.38'
    init_velocity = '${U_in} 0 0'

    gravity = ${g}
    max_substeps = 50
    max_cfl = 0.5

    debug = true

    migration_via_root = true     # default true in v4
    root_rank = 0                # router rank
    mpi_barriers = true          # optional; turn on while debugging deadlocks
    migration_tag_base = 33000   # IMPORTANT: pick a unique tag range to avoid collisions

    wall_deposition_velocity = '${U_in}'
    include_settling_in_vd = true
  []
[]

#===============================================================================
# VARIABLES
#===============================================================================
[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = ${U_in}
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    solver_sys = w_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
  []

  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = ${k_in}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = ${eps_in}
  []
[]

#===============================================================================
# FV KERNELS
#===============================================================================
[LinearFVKernels]

  # inactive = 'u_time v_time w_time TKE_time TKED_time'

  [u_time]
    type = LinearFVTimeDerivative
    variable = vel_x
    factor = ${rho}
  []
  [u_advection]
    type = LinearWCNSFVMomentumFlux
    momentum_component = 'x'
    variable = vel_x
    mu = 'mu_t'
    use_nonorthogonal_correction = true
    use_deviatoric_terms = true
  []
  [u_viscosity]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = true
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
  []

  [v_time]
    type = LinearFVTimeDerivative
    variable = vel_y
    factor = ${rho}
  []
  [v_advection]
    type = LinearWCNSFVMomentumFlux
    momentum_component = 'y'
    variable = vel_y
    mu = 'mu_t'
    use_nonorthogonal_correction = true
    use_deviatoric_terms = true
  []
  [v_viscosity]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = true
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    momentum_component = 'y'
    variable = vel_y
  []

  [w_time]
    type = LinearFVTimeDerivative
    variable = vel_z
    factor = ${rho}
  []
  [w_advection]
    type = LinearWCNSFVMomentumFlux
    momentum_component = 'z'
    variable = vel_z
    mu = 'mu_t'
    use_nonorthogonal_correction = true
    use_deviatoric_terms = true
  []
  [w_viscosity]
    type = LinearFVDiffusion
    variable = vel_z
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = true
  []
  [w_pressure]
    type = LinearFVMomentumPressure
    momentum_component = 'z'
    variable = vel_z
  []

  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = true
  []
  [p_source]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
  []

  [TKE_time]
    type = LinearFVTimeDerivative
    variable = TKE
    factor = ${rho}
  []
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = true
  []
  [TKE_diffusion_turbulent]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = true
  []
  [TKE_source_sink]
    type = LinearFVTKESourceSink
    variable = TKE
    mu = ${mu}
    epsilon = TKED
    mu_t = 'mu_t'
    walls = 'wall'
    wall_treatment = ${wall_treatment}
    # C_pl = 2.0
  []

  [TKED_time]
    type = LinearFVTimeDerivative
    variable = TKED
    factor = ${rho}
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = 'wall'
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = true
    walls = 'wall'
  []
  [TKED_diffusion_turbulent]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = true
    walls = 'wall'
  []
  [TKED_source_sink]
    type = LinearFVTKEDSourceSink
    variable = TKED
    tke = TKE
    mu = ${mu}
    mu_t = 'mu_t'
    walls = 'wall'
    wall_treatment = ${wall_treatment}
    C1_eps = ${C1_eps}
    C2_eps = 'C2_eps_functor'
    # C_pl = 2.0
  []
[]

#===============================================================================
# BOUNDARY CONDITIONS
#===============================================================================
[LinearFVBCs]
  [inlet_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_x
    functor = ${U_in}
  []
  [inlet-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_y
    functor = 0
  []
  [inlet-w]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_z
    functor = 0
  []

  # Turbulence inlet:
  # For consistency with Chen: compute TKE and TKED from I_inlet and length scale.
  # For now, we keep simple constants (tune to match Chen as needed).
  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKE
    functor = ${k_in}
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKED
    functor = ${eps_in}
  []

  [no-slip-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'wall'
    variable = vel_x
    functor = 0
  []
  [no-slip-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'wall'
    variable = vel_y
    functor = 0
  []
  [no-slip-w]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'wall'
    variable = vel_z
    functor = 0
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = 'wall'
    variable = 'mu_t'
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu  = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []

  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_x
    use_two_term_expansion = false
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_y
    use_two_term_expansion = false
  []
  [outlet_tke]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKE
    use_two_term_expansion = false
  []
  [outlet_tked]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKED
    use_two_term_expansion = false
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'outlet'
    variable = pressure
    functor = 0.0
  []
[]

#===============================================================================
# AUX VARIABLES
#===============================================================================
[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
  []

  [yplus]
    type = MooseVariableFVReal
  []

  [u_star]
    type = MooseVariableFVReal
  []

  [u_var_FV]
    type = INSFVVelocityVariable
  []
  [v_var_FV]
    type = INSFVVelocityVariable
  []
  [w_var_FV]
    type = INSFVVelocityVariable
  []

  [S2]
    type = MooseVariableFVReal
  []

  [eta]
    type = MooseVariableFVReal
  []

  [C2_eps_functor]
    type = MooseVariableFVReal
  []
[]

#===============================================================================
# AUX KERNELS
#===============================================================================
[AuxKernels]
  [compute_mu_t]
    type = kEpsilonViscosityAux
    variable = mu_t
    C_mu = ${C_mu}
    mu   = ${mu}
    epsilon = TKED
    tke = TKE
    bulk_wall_treatment = ${bulk_wall_treatment}
    walls = 'wall'
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []

  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    tke = TKE
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    walls = 'wall'
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []

  [compute_u_star]
    type = ParsedAux
    variable = u_star
    expression = 'sqrt(sqrt(${C_mu})*TKE)'
    coupled_variables = 'TKE'
    execute_on = 'NONLINEAR'
  []

  [u_var_FV]
    type = ParsedAux
    variable = u_var_FV
    expression = 'vel_x'
    coupled_variables = 'vel_x'
    execute_on = 'NONLINEAR'
  []
  [v_var_FV]
    type = ParsedAux
    variable = v_var_FV
    expression = 'vel_y'
    coupled_variables = 'vel_y'
    execute_on = 'NONLINEAR'
  []
  [w_var_FV]
    type = ParsedAux
    variable = w_var_FV
    expression = 'vel_z'
    coupled_variables = 'vel_z'
    execute_on = 'NONLINEAR'
  []

  [compute_S2]
    type = INSFVMixingLengthTurbulentViscosityAux
    variable = S2
    mixing_length = 1.0
    u = u_var_FV
    v = v_var_FV
    w = w_var_FV
    execute_on = 'NONLINEAR'
  []

  [compute_eta]
    type = ParsedAux
    variable = eta
    expression = 'sqrt(S2)*TKE/TKED'
    coupled_variables = 'S2 TKE TKED'
    execute_on = 'NONLINEAR'
  []

  [compute_C2_eps_functor]
    type = ParsedAux
    variable = C2_eps_functor
    expression = '${C2_eps} + ${C_mu}*(eta^3)*(1-eta/4.38)/(1.0 + 0.012*(eta^3))'
    coupled_variables = 'eta'
    execute_on = 'NONLINEAR'
  []
[]

#===============================================================================
# FUNCTOR MATERIALS (for drift flux & deposition BC)
#===============================================================================
[FunctorMaterials]
  # Constant fluid properties and temperature (air, isothermal)
  [fluid_properties]
    type = ADGenericFunctorMaterial
    prop_names = 'rho mu T_fluid'
    prop_values = '${rho} ${mu} 293.0'
  []
[]

#===============================================================================
# PARTICLE->FLUID COUPLING MATERIAL (force density fields)
#===============================================================================
[Materials]
  [particle_coupling]
    type = TorchParticleCouplingMaterial
    particle_cloud_userobject = particle_cloud
    output_force_density = true
    fx_name = particle_fx
    fy_name = particle_fy
    fz_name = particle_fz

    # debug = true

    outputs = exodus
    output_properties = 'aerosol_number_conc aerosol_mass_conc particle_fx particle_fy particle_fz'
  []
[]

#===============================================================================
# EXECUTIONER
#===============================================================================
[Executioner]
  type = PIMPLE

  end_time = 10.0
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    iteration_window = 2
    growth_factor = 1
    cutback_factor = 1
    dt = 0.01
  []
  num_iterations = 50
  print_fields = false

  momentum_l_max_its = 100
  continue_on_max_its = true

  momentum_systems = 'u_system v_system w_system'
  momentum_l_abs_tol = 1e-9
  momentum_l_tol     = 1e-9
  momentum_equation_relaxation = 0.65
  momentum_absolute_tolerance  = 1e-8
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'

  pressure_system = 'pressure_system'
  pressure_l_abs_tol = 1e-9
  pressure_l_tol     = 1e-9
  pressure_variable_relaxation = 0.3
  pressure_absolute_tolerance  = 1e-8
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'

  turbulence_systems = 'TKE_system TKED_system'
  turbulence_l_abs_tol = 1e-9
  turbulence_l_tol     = 1e-9
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_absolute_tolerance  = '1e-12 1e-12'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'

  # Aerosol active scalar solve group
[]

#===============================================================================
# OUTPUTS
#===============================================================================
[Outputs]
  exodus = true
  print_linear_converged_reason = false
  print_linear_residuals = false
  print_nonlinear_converged_reason = false
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
