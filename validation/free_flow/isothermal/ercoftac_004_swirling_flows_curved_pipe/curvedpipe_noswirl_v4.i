rho = 1.0
bulk_u = 10.
D = 0.0762
mu = 1.524e-5

# advected_interp_method = 'upwind'
momentum_interp_method = 'average'

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
bulk_wall_treatment = false
walls = 'wall'
wall_treatment = 'eq_newton' # Options: eq_newton, eq_incremental, eq_linearized, neq

[Mesh]
    [load_mesh]
      type = FileMeshGenerator
      #file = 'swirled_pipe_mesh_fine_cubit.e'
      file = 'curvedpipe_noswirl_lam_out.e'
      #use_for_exodus_restart = true
    []
    # [rename_blocks]
    #   type = RenameBlockGenerator
    #   input = 'load_mesh'
    #   old_block = '0 74'
    #   new_block = 'fluid fluid'
    # []
    #  [rename_boundaries]
    #    type = RenameBoundaryGenerator
    #    input = 'rename_blocks'
    #    old_boundary = '0'
    #    new_boundary = 'wall'
    #  []

    parallel_type = distributed
[]

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system TKE_system TKED_system'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    w = vel_z
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
    #initial_from_file_var = vel_x
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    #initial_condition = 0
    solver_sys = v_system
    #initial_from_file_var = vel_y
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    #initial_condition = 0
    solver_sys = w_system
    #initial_from_file_var = vel_z
  []
  [pressure]
    type = MooseLinearVariableFVReal
    #initial_condition = 1e-8
    solver_sys = pressure_system
    #initial_from_file_var = pressure
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
  [vel_z_ic]
    type = FVFunctionIC
    variable = 'vel_z'
    function = IC_vel_z
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
    advected_interp_method = ${momentum_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
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
    momentum_component = 'x'
    pressure = pressure
  []

  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    advected_interp_method = ${momentum_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
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
    momentum_component = 'y'
    pressure = pressure
  []

  [w_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_z
    advected_interp_method = ${momentum_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'z'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = true
  []
  [w_diffusion]
    type = LinearFVDiffusion
    variable = vel_z
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [w_pressure]
    type = LinearFVMomentumPressure
    variable = vel_z
    momentum_component = 'z'
    pressure = pressure
  []

  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
  []
  [p_source]
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
  [TKE_diffusion_turbulent]
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
    w = vel_z
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
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
  [TKED_diffusion_turbulent]
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
    w = vel_z
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
  [inlet-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_x
    functor = 0
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
    #functor = ${bulk_u}
    functor = 'fully_developed_velocity'
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

  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'outlet'
    variable = pressure
    functor = 0
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
  [outlet_w]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = vel_z
    use_two_term_expansion = false
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKE
    use_two_term_expansion = false
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKED
    use_two_term_expansion = false
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
  [walls-w]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_z
    functor = 0
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = ${walls}
    variable = mu_t
    u = vel_x
    v = vel_y
    w = vel_z
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
    prop_file_name = 'FDFlow_test.csv'
    read_type = 'voronoi'
    #nprop = 13 # number of columns in CSV
    nprop = 6
    execute_on = TIMESTEP_BEGIN
    nvoronoi = 1125
  []
  [read_IC]
    type = PropertyReadFile
    prop_file_name = 'IC_3000eq.csv'
    read_type = 'voronoi'
    nprop = 9 # number of columns in CSV
    execute_on = TIMESTEP_BEGIN
    nvoronoi = 540000
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
    column_number = '4'
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
  [IC_vel_z]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = 'read_IC'
    read_type = 'voronoi'
    column_number = '8'
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
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
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
    w = vel_z
    walls = ${walls}
    bulk_wall_treatment = ${bulk_wall_treatment}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e20
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
    w = vel_z
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []
[]


[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.25 0.25'
  turbulence_field_relaxation = '0.25 0.25'
  num_iterations = 1500
  pressure_absolute_tolerance = 1e-12
  momentum_absolute_tolerance = 1e-12
  turbulence_absolute_tolerance = '1e-12 1e-12'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'

  momentum_l_abs_tol = 1e-14
  pressure_l_abs_tol = 1e-14
  turbulence_l_abs_tol = 1e-14
  momentum_l_max_its = 30
  pressure_l_max_its = 30
  turbulence_l_max_its = 30
  momentum_l_tol = 0.0
  pressure_l_tol = 0.0
  turbulence_l_tol = 0.0
  print_fields = false

  continue_on_max_its = true
[]

[Outputs]
  [out]
    type = Exodus
  []
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
