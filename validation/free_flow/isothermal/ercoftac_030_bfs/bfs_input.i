########################################################################
# BFS simulation
########################################################################

# Modeling Parameters
rho = 1.18415
bulk_u = 44.2
D = 0.1016
mu = 1.8551e-5

advected_interp_method = 'upwind'

### k-epsilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09
walls = 'wall'
wall_treatment = 'neq' # Options: eq_newton, eq_incremental, eq_linearized, neq

### Initial and Boundary Conditions ###
intensity = 0.01
k_init = '${fparse 1.5*(intensity * bulk_u)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / D}'

### Postprocessing
y_first_cell_1 = 5e-4
y_first_cell_2 = ${fparse (0.0127 - 0.0119) / 2}

# Turbulence-model knobs (optional)
k_epsilon_variant   = 'Standard'  # Standard | StandardLowRe | StandardTwoLayer | Realizable | RealizableTwoLayer
# two_layer_flavor    = 'Wolfstein' # Wolfstein | NorrisReynolds | Xu (only used for *TwoLayer variants)
# use_buoyancy        = false
# use_compressibility = false
nonlinear_model     = 'none'
# curvature_model     = 'none'
use_yap             = false
# use_low_re_Gprime   = false
bulk_wall_treatment = false

[Mesh]
  [load_mesh]
    type = FileMeshGenerator
    file = 'BFS_mesh_fine1.e'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    pressure_projection_method = 'consistent'
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
  []
[]

[FVICs]
  [vel_x_ic]
     type = FVFunctionIC
     variable = 'vel_x'
     function = fully_developed_velocity
  []
  [vel_y_ic]
    type = FVFunctionIC
    variable = 'vel_y'
    function = 0
  []
  [pressure_ic]
    type = FVFunctionIC
    variable = 'pressure'
    function = 0
  []
  [TKE_ic]
    type = FVFunctionIC
    variable = 'TKE'
    function = fully_developed_tke
  []
  [TKED_ic]
    type = FVFunctionIC
    variable = 'TKED'
    function = fully_developed_tked
  []
[]

[LinearFVKernels]
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = true
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = true
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
  []
  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
  []
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink]
    type = kEpsilonTKESourceSink
    variable = TKE

    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 2

    # NEW (optional)
    k_epsilon_variant        = ${k_epsilon_variant}    # 'Standard', 'Realizable', etc.
    nonlinear_model          = ${nonlinear_model}

    # use_buoyancy             = ${use_buoyancy}
    # use_compressibility      = ${use_compressibility}
    # curvature_model          = ${curvature_model}
    # Pr_t                     = 0.9
    # C_M                      = 1.0
  []
  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_source_sink]
    type = kEpsilonTKEDSourceSink
    variable = TKED

    u = vel_x
    v = vel_y
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 2

    # NEW (optional)
    k_epsilon_variant   = ${k_epsilon_variant}
    use_yap             = ${use_yap}
    nonlinear_model     = ${nonlinear_model}
    wall_distance       = distance

    # use_buoyancy        = ${use_buoyancy}
    # use_compressibility = ${use_compressibility}
    # curvature_model     = ${curvature_model}
    # use_low_re_Gprime   = ${use_low_re_Gprime}
    # Pr_t = 0.9
    # C_M  = 1.0
  []
[]

[LinearFVBCs]
  [inlet_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_x
    functor = 'fully_developed_velocity'
  []
  [inlet-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_y
    functor = 0
  []
  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKE
    functor = 'fully_developed_tke'
  []
  # Alternative constant inlet BC
  # [inlet_TKE]
  #   type = INSFVInletIntensityTKEBC
  #   boundary = 'inlet'
  #   variable = TKE
  #   u = vel_x
  #   v = vel_y
  #   w = vel_z
  #   intensity = ${intensity}
  # []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKED
    functor = 'fully_developed_tked'
  []
  # Alternative constant inlet BC
  # [inlet_TKED]
  #   type = INSFVMixingLengthTKEDBC
  #   boundary = 'inlet'
  #   variable = TKED
  #   k = TKE
  #   characteristic_length = '${D}'
  # []
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
  [walls-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_x
    functor = 0
  []
  [walls-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_y
    functor = 0
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = ${walls}
    variable = 'mu_t'
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
[]

[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
  []
  # For computing c_f
  [mu_t_wall]
    type = MooseLinearVariableFVReal
  []
  [yplus]
    type = MooseLinearVariableFVReal
  []
  [mu_eff]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse mu + rho * C_mu * ${k_init}^2 / eps_init}'
  []
  [distance]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [compute_mu_t]
    type = kEpsilonViscosity
    variable = mu_t

    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y

    bulk_wall_treatment = ${bulk_wall_treatment}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e20
    execute_on = 'NONLINEAR'

    # NEW (optional) – choose model and options
    k_epsilon_variant = ${k_epsilon_variant}    # e.g. 'Standard' or 'Realizable'
    # two_layer_flavor  = ${two_layer_flavor}     # ignored unless *TwoLayer variants
    # Cd0 = 0.091      # defaults, can omit if you keep Standard
    # Cd1 = 0.0042
    # Cd2 = 0.00011
    # Ca0 = 0.667      # Realizable C_mu coefficients
    # Ca1 = 1.25
    # Ca2 = 1.0
    # Ca3 = 0.9
    wall_distance = distance   # only needed for LowRe/TwoLayer (see below)
  []
  [compute_mu_t_wall_neq]
    type = ParsedAux
    variable = 'mu_t_wall'
    expression = '${mu} * (0.4187 * yplus / log(max(9.793 * yplus, 1.0 + 1e-4)) - 1)'
    coupled_variables = 'yplus'
    execute_on = 'TIMESTEP_END'
  []
  # Equivalent way of computing mu_t_wall
  # [compute_mu_t_wall]
  #   type = kEpsilonViscosityAux
  #   variable = mu_t_wall
  #   C_mu = ${C_mu}
  #   tke = TKE
  #   epsilon = TKED
  #   mu = ${mu}
  #   rho = ${rho}
  #   u = vel_x
  #   v = vel_y
  #   bulk_wall_treatment = true
  #   walls = ${walls}
  #   wall_treatment = ${wall_treatment}
  #   execute_on = 'TIMESTEP_END'
  #   mu_t_ratio_max = 1e20
  # []
  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    tke = TKE
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []
  [compute_mu_eff]
    type = ParsedAux
    variable = 'mu_eff'
    coupled_variables = 'mu_t'
    expression = 'mu_t + ${mu}'
    execute_on = 'NONLINEAR'
  []
  [distance_aux]
    type = WallDistanceMixingLengthAux
    walls = 'wall'
    variable = 'distance'
    execute_on = 'initial'
    von_karman_const = 1.0
  []
[]

[UserObjects]
  [read_recycling]
    type = PropertyReadFile
    prop_file_name = 'BFS_FDprofile.csv'
    read_type = 'voronoi'
    nprop = 7
    execute_on = TIMESTEP_BEGIN
    nvoronoi = 60
  []
[]

[Functions]
  [fully_developed_velocity]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_recycling'
    read_type = 'voronoi'
    column_number = '6'
  []
  [fully_developed_tke]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_recycling'
    read_type = 'voronoi'
    column_number = '3'
  []
  [fully_developed_tked]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_recycling'
    read_type = 'voronoi'
    column_number = '4'
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  momentum_l_abs_tol = 1e-8
  pressure_l_abs_tol = 1e-8
  turbulence_l_abs_tol = 1e-8
  momentum_l_tol = 1e-10
  pressure_l_tol = 1e-10
  turbulence_l_tol = 1e-10
  momentum_equation_relaxation = 0.9
  pressure_variable_relaxation = 0.5
  turbulence_equation_relaxation = '0.2 0.2'
  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  turbulence_absolute_tolerance = '1e-8 1e-8'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  momentum_l_max_its = 300
  pressure_l_max_its = 300
  turbulence_l_max_its = 30
  print_fields = false

  num_iterations = 6000
  continue_on_max_its = true
[]

[VectorPostprocessors]
  [inlet_channel_wall_sampler]
    type = LineValueSampler
    variable = 'distance vel_x pressure mu_t_wall'
    start_point = '${fparse -0.048} ${y_first_cell_1} 0.0'
    end_point   = '${fparse -1e-3}  ${y_first_cell_1} 0.0'
    sort_by = 'x'
    num_points = 100
  []
  [outlet_channel_wall_sampler]
    type = LineValueSampler
    variable = 'distance vel_x pressure mu_t_wall'
    start_point = '6.1e-4 ${fparse -0.0127 + y_first_cell_2} 0.0'
    end_point =   '0.2492 ${fparse -0.0127 + y_first_cell_2} 0.0'
    sort_by = 'x'
    num_points = 100
  []
  [vel_x_xoH_1_sampler]
    type = LineValueSampler
    variable = 'vel_x'
    start_point = '${fparse 0.0127*1.0} -0.0127 0.0'
    end_point = '${fparse 0.0127*1.0} ${fparse 0.0127*8.2-0.0127} 0.0'
    sort_by = 'y'
    num_points = 100
  []
  [vel_x_xoH_4_sampler]
    type = LineValueSampler
    variable = 'vel_x'
    start_point = '${fparse 0.0127*4.0} -0.0127 0.0'
    end_point = '${fparse 0.0127*4.0} ${fparse 0.0127*8.2-0.0127} 0.0'
    sort_by = 'y'
    num_points = 100
  []
  [vel_x_xoH_6_sampler]
    type = LineValueSampler
    variable = 'vel_x'
    start_point = '${fparse 0.0127*6.0} -0.0127 0.0'
    end_point = '${fparse 0.0127*6.0} ${fparse 0.0127*8.2-0.0127} 0.0'
    sort_by = 'y'
    num_points = 100
  []
  [vel_x_xoH_10_sampler]
    type = LineValueSampler
    variable = 'vel_x'
    start_point = '${fparse 0.0127*10.0} -0.0127 0.0'
    end_point = '${fparse 0.0127*10.0} ${fparse 0.0127*8.2-0.0127} 0.0'
    sort_by = 'y'
    num_points = 100
  []
[]

[Postprocessors]
  [mdot]
    type = VolumetricFlowRate
    vel_x = vel_x
    vel_y = vel_y
    advected_interp_method = ${advected_interp_method}
    advected_quantity = ${rho}
    boundary = 'inlet'
    execute_on = 'TIMESTEP_END'
  []
  [bulk_u]
    type = SideExtremeValue
    boundary = 'inlet'
    variable = 'vel_x'
    execute_on = 'TIMESTEP_END'
  []
  [Reynolds_in]
    type = ParsedPostprocessor
    expression = 'bulk_u / mu * rho * D'
    pp_names = 'bulk_u'
    constant_expressions = '${mu} ${rho} ${D}'
    constant_names = 'mu rho D'
    execute_on = 'TIMESTEP_END'
  []
[]

[Outputs]
  [exo]
    type = Exodus
  []
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
