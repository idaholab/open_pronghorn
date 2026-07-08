rho = 1.0
bulk_u = 10.
D = 0.0762
mu = 1.524e-5

advected_interp_method = 'upwind'
momentum_interp_method = 'average'

pressure_tag = "pressure_grad"

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
[]

[Problem]
  nl_sys_names = 'u_system v_system w_system pressure_system TKE_system TKED_system'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  #advected_interp_method = ${advected_interp_method}
  velocity_interp_method = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolatorSegregated
    u = vel_x
    v = vel_y
    w = vel_z
    pressure = pressure
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    solver_sys = u_system
    two_term_boundary_expansion = false
    #initial_from_file_var = vel_x
  []
  [vel_y]
    type = INSFVVelocityVariable
    #initial_condition = 0
    solver_sys = v_system
    two_term_boundary_expansion = false
    #initial_from_file_var = vel_y
  []
  [vel_z]
    type = INSFVVelocityVariable
    #initial_condition = 0
    solver_sys = w_system
    two_term_boundary_expansion = false
    #initial_from_file_var = vel_z
  []
  [pressure]
    type = INSFVPressureVariable
    #initial_condition = 1e-8
    solver_sys = pressure_system
    two_term_boundary_expansion = false
    #initial_from_file_var = pressure
  []
  [TKE]
    type = INSFVEnergyVariable
    solver_sys = TKE_system
    #initial_condition = ${k_init}
  []
  [TKED]
    type = INSFVEnergyVariable
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

[FVKernels]
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
    advected_interp_method = ${momentum_interp_method}
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_viscosity_turbulent]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = 'mu_t'
    momentum_component = 'x'
    complete_expansion = true
    u = vel_x
    v = vel_y
    w = vel_z
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
    extra_vector_tags = ${pressure_tag}
  []

  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
    advected_interp_method = ${momentum_interp_method}
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_viscosity_turbulent]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = 'mu_t'
    momentum_component = 'y'
    complete_expansion = true
    u = vel_x
    v = vel_y
    w = vel_z
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
    extra_vector_tags = ${pressure_tag}
  []

  [w_advection]
    type = INSFVMomentumAdvection
    variable = vel_z
    rho = ${rho}
    momentum_component = 'z'
    advected_interp_method = ${momentum_interp_method}
  []
  [w_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_z
    mu = ${mu}
    momentum_component = 'z'
  []
  [w_viscosity_turbulent]
    type = INSFVMomentumDiffusion
    variable = vel_z
    mu = 'mu_t'
    momentum_component = 'z'
    complete_expansion = true
    u = vel_x
    v = vel_y
    w = vel_z
  []
  [w_pressure]
    type = INSFVMomentumPressure
    variable = vel_z
    momentum_component = 'z'
    pressure = pressure
    extra_vector_tags = ${pressure_tag}
  []

  [p_diffusion]
    type = FVAnisotropicDiffusion
    variable = pressure
    coeff = "Ainv"
    coeff_interp_method = 'average'
  []
  [p_source]
    type = FVDivergence
    variable = pressure
    vector_field = "HbyA"
    force_boundary_execution = true
  []

  [TKE_advection]
    type = INSFVTurbulentAdvection
    variable = TKE
    rho = ${rho}
    advected_interp_method = ${momentum_interp_method}
  []
  [TKE_diffusion]
    type = INSFVTurbulentDiffusion
    variable = TKE
    coeff = ${mu}
  []
  [TKE_diffusion_turbulent]
    type = INSFVTurbulentDiffusion
    variable = TKE
    coeff = 'mu_t'
    scaling_coef = ${sigma_k}
  []
  [TKE_source_sink]
    type = INSFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    w = vel_z
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
  []

  [TKED_advection]
    type = INSFVTurbulentAdvection
    variable = TKED
    rho = ${rho}
    walls = ${walls}
    advected_interp_method = ${momentum_interp_method}
  []
  [TKED_diffusion]
    type = INSFVTurbulentDiffusion
    variable = TKED
    coeff = ${mu}
    walls = ${walls}
  []
  [TKED_diffusion_turbulent]
    type = INSFVTurbulentDiffusion
    variable = TKED
    coeff = 'mu_t'
    scaling_coef = ${sigma_eps}
    walls = ${walls}
   []
  [TKED_source_sink]
    type = INSFVTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    w = vel_z
    k = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = vel_x
    functor = 0
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = vel_y
    functor = 0
  []
  [inlet-w]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = vel_z
    #functor = ${bulk_u}
    functor = 'fully_developed_velocity'
  []
  [inlet_TKE]
    type = FVFunctionDirichletBC
    boundary = 'inlet'
    variable = TKE
    function = 'fully_developed_tke'
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
    type = FVFunctionDirichletBC
    boundary = 'inlet'
    variable = TKED
    function = 'fully_developed_tked'
  []
  # [inlet_TKED]
  #   type = INSFVMixingLengthTKEDBC
  #   boundary = 'inlet'
  #   variable = TKED
  #   k = TKE
  #   characteristic_length = '${D}'
  # []

  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'outlet'
    variable = pressure
    functor = 0
  []

  [walls-u]
    type = FVDirichletBC
    boundary = ${walls}
    variable = vel_x
    value = 0
  []
  [walls-v]
    type = FVDirichletBC
    boundary = ${walls}
    variable = vel_y
    value = 0
  []
  [walls-w]
    type = FVDirichletBC
    boundary = ${walls}
    variable = vel_z
    value = 0
  []
  [walls_mu_t]
    type = INSFVTurbulentViscosityWallFunction
    boundary = ${walls}
    variable = mu_t
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    k = TKE
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
    type = MooseVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
    two_term_boundary_expansion = false
  []
  [yplus]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
  []
[]

[AuxKernels]
  [compute_mu_t]
    type = kEpsilonViscosityAux
    variable = mu_t
    C_mu = ${C_mu}
    k = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    w = vel_z
    walls = ${walls}
    execute_on = 'NONLINEAR'
  []
  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    k = TKE
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
  type = SIMPLENonlinearAssembly
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKED_system TKE_system'

  pressure_gradient_tag = ${pressure_tag}
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.25 0.25'
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
  momentum_l_tol = 0.0
  pressure_l_tol = 0.0
  turbulence_l_tol = 0.0
  print_fields = false
[]

[Outputs]
  [v4_out]
    type = Exodus
  []
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]
