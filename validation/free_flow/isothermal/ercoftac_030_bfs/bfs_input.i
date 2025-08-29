
# Modeling Parameters
rho = 1.18415
bulk_u = 44.2
D = 0.2032
mu = 1.8551e-5

advected_interp_method = 'upwind'

### k-epslilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09

### Initial and Boundary Conditions ###
intensity = 0.01
k_init = '${fparse 1.5*(intensity * bulk_u)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / D}'
### Modeling parameters ###
walls = 'wall'
wall_treatment = 'neq' # Options: eq_newton, eq_incremental, eq_linearized, neq
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
  #velocity_interp_method = 'rc'
[]
[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
  []
[]
[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    #initial_condition = ${bulk_u}
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    #initial_condition = 0
    solver_sys = v_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    #initial_condition = 1e-8
    solver_sys = pressure_system
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    #initial_condition = ${k_init}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    #initial_condition = ${eps_init}
  []
[]
[FVICs]
  [vel_x_ic]
     type = FVFunctionIC
     variable = 'vel_x'
     function = IC_vel_x
  []
  [vel_y_ic]
    type = FVFunctionIC
    variable = 'vel_y'
    function = IC_vel_y
  []
  [pressure_ic]
    type = FVFunctionIC
    variable = 'pressure'
    function = IC_pressure
  []
  [TKE_ic]
    type = FVFunctionIC
    variable = 'TKE'
    function = IC_TKE
  []
  [TKED_ic]
    type = FVFunctionIC
    variable = 'TKED'
    function = IC_TKED
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
    use_nonorthogonal_correction = true
    use_deviatoric_terms = false
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = true
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
    use_nonorthogonal_correction = true
    use_deviatoric_terms = false
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = true
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
    use_nonorthogonal_correction = true
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
    use_nonorthogonal_correction = true
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = true
  []
  [TKE_source_sink]
    type = LinearFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 2.0
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
    use_nonorthogonal_correction = true
    walls = ${walls}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = true
    walls = ${walls}
  []
  [TKED_source_sink]
    type = LinearFVTKEDSourceSink
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
[UserObjects]
  [read_recycling]
    type = PropertyReadFile
    #prop_file_name = 'FDflow.csv'
    prop_file_name = 'BFS_FDprofile.csv'
    read_type = 'voronoi'
    #nprop = 13 # number of columns in CSV
    nprop = 7
    execute_on = TIMESTEP_BEGIN
    nvoronoi = 60
  []
  [read_IC]
    type = PropertyReadFile
    prop_file_name = 'BFS_IC.csv'
    read_type = 'voronoi'
    nprop = 8 # number of columns in CSV
    execute_on = TIMESTEP_BEGIN
    nvoronoi = 10575
  []
[]
[Functions]
  [fully_developed_velocity]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_recycling'
    read_type = 'voronoi'
    column_number = '3'
  []
  [fully_developed_tke]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_recycling'
    read_type = 'voronoi'
    column_number = '6'
  []
  [fully_developed_tked]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_recycling'
    read_type = 'voronoi'
    column_number = '5'
  []
  [IC_vel_x]
   type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_IC'
    read_type = 'voronoi'
    column_number = '6'
  []
  [IC_vel_y]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_IC'
    read_type = 'voronoi'
    column_number = '7'
  []
  [IC_pressure]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_IC'
    read_type = 'voronoi'
    column_number = '3'
  []
  [IC_TKE]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_IC'
    read_type = 'voronoi'
    column_number = '4'
  []
  [IC_TKED]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_IC'
    read_type = 'voronoi'
    column_number = '5'
  []
[]
[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
  []
  [yplus]
    type = MooseLinearVariableFVReal
  []
  [mu_eff]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
  []
  [distance]
    type = MooseVariableFVReal
  []
[]
[AuxKernels]
  [compute_mu_t]
    type = kEpsilonViscosityAux
    variable = mu_t
    C_mu = ${C_mu}
    tke =TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    bulk_wall_treatment = false
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
    mu_t_ratio_max = 1e20
  []
  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    tke =TKE
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
[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  momentum_l_abs_tol = 1e-14
  pressure_l_abs_tol = 1e-14
  turbulence_l_abs_tol = 1e-14
  momentum_l_tol = 1e-14
  pressure_l_tol = 1e-14
  turbulence_l_tol = 1e-14
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.25 0.25'
  num_iterations = 3000
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
  continue_on_max_its = true
[]
[VectorPostprocessors]
  [inlet_sampler]
    type = LineValueSampler
    variable = 'distance vel_x pressure mu_t'
    start_point = '-0.05 1e-5 0.0'
    end_point = '1e-5 1e-5 0.0'
    sort_by = 'x'
    num_points = 100
  []
  [outlet_sampler]
    type = LineValueSampler
    variable = 'distance vel_x pressure mu_t'
    start_point = '0.0 ${fparse -0.0127 + 1e-5} 0.0'
    end_point = '0.25 ${fparse -0.0127 + 1e-5} 0.0'
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
[Outputs]
  [durbin]
    type = Exodus
  []
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
