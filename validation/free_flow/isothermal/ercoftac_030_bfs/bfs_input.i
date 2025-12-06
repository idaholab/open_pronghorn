########################################################################
# BFS simulation
########################################################################

# Modeling Parameters
rho = 1.18415
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
bulk_wall_treatment = false

### Initial Conditions for Non-Developed Case ###
# intensity = 0.01
# bulk_u = 44.2
# k_init = '${fparse 1.5*(intensity * bulk_u)^2}'
# eps_init = '${fparse C_mu^0.75 * k_init^1.5 / D}'

### Postprocessing
y_first_cell_1 = 5e-4
y_first_cell_2 = ${fparse (0.0127 - 0.0119) / 2}

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

[Physics]
  [NavierStokes]
    [FlowSegregated/flow]
      velocity_variable = 'vel_x vel_y'
      pressure_variable = 'pressure'

      # Initial conditions
      initial_velocity = 'fully_developed_velocity 0 0'
      initial_pressure = '1e-8'

      # Material properties
      density = ${rho}
      dynamic_viscosity = ${mu}

      # Boundary conditions
      inlet_boundaries = 'inlet'
      momentum_inlet_types = 'fixed-velocity'
      momentum_inlet_functors = 'fully_developed_velocity 0'

      wall_boundaries = ${walls}
      momentum_wall_types = 'noslip'

      outlet_boundaries = 'outlet'
      momentum_outlet_types = 'fixed-pressure'
      pressure_functors = '0'

      # Numerical parameters
      include_deviatoric_stress = true
      orthogonality_correction = false
      pressure_two_term_bc_expansion = false
      momentum_two_term_bc_expansion = false
      momentum_advection_interpolation = ${advected_interp_method}
      system_names = 'u_system v_system pressure_system'
      # consistency with kernel syntax input, allows to use the same executioner block
      rhie_chow_uo_name = 'rc'
      pressure_projection_method = 'consistent'
    []
    [TurbulenceSegregated/k-epsilon]
      # Model
      turbulence_handling = 'k-epsilon'
      tke_name = TKE
      tked_name = TKED
      system_names = 'TKE_system TKED_system'

      initial_tke = 'fully_developed_tke'
      initial_tked = 'fully_developed_tked'
      initial_mu_t = 'initial_mu_t'

      # Model parameters
      mu_t_ratio_max = 1e20
      sigma_k = ${sigma_k}
      sigma_eps = ${sigma_eps}
      C_pl = 2
      C1_eps = ${C1_eps}
      C2_eps = ${C2_eps}

      turbulence_walls = ${walls}
      wall_treatment_eps = ${wall_treatment}
      bulk_wall_treatment = ${bulk_wall_treatment}
      use_nonorthogonal_correction = false
    []
  []
[]

[LinearFVBCs]
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
[]

[AuxVariables]
  # For computing c_f
  [mu_t_wall]
    type = MooseLinearVariableFVReal
  []
  [yplus]
    type = MooseLinearVariableFVReal
  []
  [distance]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
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
  [initial_mu_t]
    type = ParsedFunction
    expression = '${rho} * ${C_mu} * TKE_initial^2 / TKED_initial'
    symbol_names = 'TKE_initial TKED_initial'
    symbol_values = 'fully_developed_tke fully_developed_tked'
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
    rhie_chow_user_object = 'rc'
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
