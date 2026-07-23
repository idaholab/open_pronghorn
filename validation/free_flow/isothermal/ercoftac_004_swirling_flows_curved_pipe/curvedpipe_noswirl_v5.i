rho = 1.0
bulk_u = 10.4 #10.
D = 0.0762
mu = 1.58e-5 #1.524e-5

# advected_interp_method = 'upwind'
momentum_interp_method = 'upwind'     

### k-epslilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09
C_pl = 10

### Initial and Boundary Conditions ###
intensity = 0.01
k_init = '${fparse 1.5*(intensity * bulk_u)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / D}'

### Modeling parameters ###
bulk_wall_treatment = false
walls = 'wall'
wall_treatment = 'neq' # Options: eq_newton, eq_incremental (best), eq_linearized, neq
k_epsilon_variant   = 'Realizable'  # Standard | StandardLowRe | StandardTwoLayer | Realizable | RealizableTwoLayer
# two_layer_flavor    = 'Wolfstein' # Wolfstein | NorrisReynolds | Xu (only used for *TwoLayer variants)
use_buoyancy        = false
use_compressibility = false
nonlinear_model     = 'none'
curvature_model     = 'standard'
use_yap             = false  #to try later on
use_low_re_Gprime   = false
use_curvature_correction = true

[Mesh]
    [load_mesh]
      type = FileMeshGenerator
      #file = 'swirled_pipe_mesh_fine_cubit.e'
      # file = 'curvedpipe_noswirl_lam_out.e'
      file = 'curvedpipe_noswirl_lam_wall_coarsened.e'
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
    C_pl = ${C_pl}
    k_epsilon_variant        = ${k_epsilon_variant}
    use_buoyancy             = ${use_buoyancy}
    use_compressibility      = ${use_compressibility}
    nonlinear_model          = ${nonlinear_model}
    curvature_model          = ${curvature_model}
    use_curvature_correction = ${use_curvature_correction}

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
    C_pl = ${C_pl}
    k_epsilon_variant   = ${k_epsilon_variant}
    use_buoyancy        = ${use_buoyancy}
    use_compressibility = ${use_compressibility}
    nonlinear_model     = ${nonlinear_model}
    curvature_model     = ${curvature_model}
    use_yap             = ${use_yap}
    use_low_re_Gprime   = ${use_low_re_Gprime}
    wall_distance       = distance
    use_curvature_correction = ${use_curvature_correction}
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
    # prop_file_name = 'FDFlow_sm18.csv'
    # prop_file_name = 'FDFlow_sm18_reconstructed.csv'
    prop_file_name = 'FDFlow_sm02_reconstructed.csv'
    # prop_file_name = 'FDFlow_test.csv'
    # prop_file_name = 'FDFlow_test_sm02interp.csv'
    # prop_file_name = 'FDFlow_test_sm02interp_massnorm.csv'
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
  [mu_t_wall]
    type = MooseLinearVariableFVReal
  []
  [yplus]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
  []
  [mu_eff]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse mu + rho * C_mu * ${k_init}^2 / eps_init}'
  []
  # Standard k-epsilon does not use this, but the BFS-style model hooks accept it.
  # Do not compute it here: the available wall-distance aux kernels require replicated meshes.
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
    w = vel_z
    walls = ${walls}
    bulk_wall_treatment = ${bulk_wall_treatment}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e20
    execute_on = 'NONLINEAR'
    k_epsilon_variant = ${k_epsilon_variant}
    wall_distance = distance
  []
  [compute_mu_t_wall_neq]
    type = ParsedAux
    variable = 'mu_t_wall'
    expression = '${mu} * (0.4187 * yplus / log(max(9.793 * yplus, 1.0 + 1e-4)) - 1)'
    coupled_variables = 'yplus'
    execute_on = 'TIMESTEP_END'
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
  [compute_mu_eff]
    type = ParsedAux
    variable = 'mu_eff'
    coupled_variables = 'mu_t'
    expression = 'mu_t + ${mu}'
    execute_on = 'NONLINEAR'
  []
[]


[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'

  momentum_equation_relaxation = 0.4
  pressure_variable_relaxation = 0.1
  turbulence_equation_relaxation = '0.25 0.25'
  turbulence_field_relaxation = '0.25 0.25'
  num_iterations = 5000 # 1000
  pressure_absolute_tolerance = 1e-15
  momentum_absolute_tolerance = 1e-15
  turbulence_absolute_tolerance = '1e-15 1e-15'
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

### Sampling stations ###
# Geometry derived directly from curvedpipe_noswirl_lam_out.e (verified via circle fit
# on the wall node coordinates, not assumed):
#   pipe radius a = 0.0381 m (D = 0.0762 m)
#   bend center   = (x=-0.4953, y=0, z=0), bend radius R = 0.4953 m, bend plane = X-Z
#   outlet leg    = straight run at x=-0.9906, z in [-1.3716, 0], flow in -z direction
#   => s/D = n downstream of the bend exit (z=0) is at z = -n*D
#   y is the "vertical plane BB" direction used for the paper's Fig. 9 comparison
[VectorPostprocessors]
  [bev00_sp01]
    type = PointValueSampler
    points = '-0.9906 0.0000000 -0.0762
              -0.9906 0.0047625 -0.0762
              -0.9906 0.0095250 -0.0762
              -0.9906 0.0142875 -0.0762
              -0.9906 0.0190500 -0.0762
              -0.9906 0.0238125 -0.0762
              -0.9906 0.0269748 -0.0762
              -0.9906 0.0303657 -0.0762
              -0.9906 0.0333375 -0.0762
              -0.9906 0.0357378 -0.0762'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []
  
  [bev00_sm01]
    type = PointValueSampler
    points = '0.0 0.0000000 -0.0762
              0.0 0.0047625 -0.0762
              0.0 0.0095250 -0.0762
              0.0 0.0142875 -0.0762
              0.0 0.0190500 -0.0762
              0.0 0.0238125 -0.0762
              0.0 0.0269748 -0.0762
              0.0 0.0303657 -0.0762
              0.0 0.0333375 -0.0762
              0.0 0.0357378 -0.0762'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  # --- Bend stations (vertical plane BB, y-direction, +y side only) ---
  # Centerline at each theta: (Cx + R cos(theta), 0, Cz + R sin(theta))
  # with Cx=-0.4953, Cz=0, R=0.4953 (verified against wall mesh, ~0.002-0.003 m match)
  [bev00_th0225]
    type = PointValueSampler
    points = '-0.0377025 0.0000000 0.1895431
              -0.0377025 0.0047625 0.1895431
              -0.0377025 0.0095250 0.1895431
              -0.0377025 0.0142875 0.1895431
              -0.0377025 0.0190500 0.1895431
              -0.0377025 0.0238125 0.1895431
              -0.0377025 0.0269748 0.1895431
              -0.0377025 0.0303657 0.1895431
              -0.0377025 0.0333375 0.1895431
              -0.0377025 0.0357378 0.1895431'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  [bev00_th0675]
    type = PointValueSampler
    points = '-0.3057569 0.0000000 0.4575975
              -0.3057569 0.0047625 0.4575975
              -0.3057569 0.0095250 0.4575975
              -0.3057569 0.0142875 0.4575975
              -0.3057569 0.0190500 0.4575975
              -0.3057569 0.0238125 0.4575975
              -0.3057569 0.0269748 0.4575975
              -0.3057569 0.0303657 0.4575975
              -0.3057569 0.0333375 0.4575975
              -0.3057569 0.0357378 0.4575975'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  [bev00_th1125]
    type = PointValueSampler
    points = '-0.6848431 0.0000000 0.4575975
              -0.6848431 0.0047625 0.4575975
              -0.6848431 0.0095250 0.4575975
              -0.6848431 0.0142875 0.4575975
              -0.6848431 0.0190500 0.4575975
              -0.6848431 0.0238125 0.4575975
              -0.6848431 0.0269748 0.4575975
              -0.6848431 0.0303657 0.4575975
              -0.6848431 0.0333375 0.4575975
              -0.6848431 0.0357378 0.4575975'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  [bev00_th1575]
    type = PointValueSampler
    points = '-0.9528975 0.0000000 0.1895431
              -0.9528975 0.0047625 0.1895431
              -0.9528975 0.0095250 0.1895431
              -0.9528975 0.0142875 0.1895431
              -0.9528975 0.0190500 0.1895431
              -0.9528975 0.0238125 0.1895431
              -0.9528975 0.0269748 0.1895431
              -0.9528975 0.0303657 0.1895431
              -0.9528975 0.0333375 0.1895431
              -0.9528975 0.0357378 0.1895431'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  # --- Downstream leg stations (vertical plane BB), outlet leg x=-0.9906 ---
  [bev00_sp06]
    type = PointValueSampler
    points = '-0.9906000 0.0000000 -0.4572000
              -0.9906000 0.0047625 -0.4572000
              -0.9906000 0.0095250 -0.4572000
              -0.9906000 0.0142875 -0.4572000
              -0.9906000 0.0190500 -0.4572000
              -0.9906000 0.0238125 -0.4572000
              -0.9906000 0.0269748 -0.4572000
              -0.9906000 0.0303657 -0.4572000
              -0.9906000 0.0333375 -0.4572000
              -0.9906000 0.0357378 -0.4572000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  [bev00_sp10]
    type = PointValueSampler
    points = '-0.9906000 0.0000000 -0.7620000
              -0.9906000 0.0047625 -0.7620000
              -0.9906000 0.0095250 -0.7620000
              -0.9906000 0.0142875 -0.7620000
              -0.9906000 0.0190500 -0.7620000
              -0.9906000 0.0238125 -0.7620000
              -0.9906000 0.0269748 -0.7620000
              -0.9906000 0.0303657 -0.7620000
              -0.9906000 0.0333375 -0.7620000
              -0.9906000 0.0357378 -0.7620000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  # NOTE: this station sits exactly at the outlet boundary plane (z=-1.3716,
  # the domain's outflow BC location) -- not a fully independent interior
  # check, see caveat in chat.
  [bev00_sp18]
    type = PointValueSampler
    points = '-0.9906000 0.0000000 -1.3716000
              -0.9906000 0.0047625 -1.3716000
              -0.9906000 0.0095250 -1.3716000
              -0.9906000 0.0142875 -1.3716000
              -0.9906000 0.0190500 -1.3716000
              -0.9906000 0.0238125 -1.3716000
              -0.9906000 0.0269748 -1.3716000
              -0.9906000 0.0303657 -1.3716000
              -0.9906000 0.0333375 -1.3716000
              -0.9906000 0.0357378 -1.3716000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'y'
    execute_on = 'FINAL'
  []

  # --- Horizontal plane (AA, in bend-plane, FULL DIAMETER) sampler blocks ---
  # Sign convention (r/a as in the ERCOFTAC .dat files) verified empirically against
  # the Dean-effect axial-velocity asymmetry in beh00-sp01.dat / beh00-th1575.dat:
  # NEGATIVE r/a = outer bend wall (high W), POSITIVE r/a = inner bend wall (low W).
  # Point formula: Global = C + (R - r_dat)*(cos(theta), 0, sin(theta)) on the bend,
  # with straight legs as the theta=0 / theta=180 limiting case plus axial translation.
  # All 19 points per station (full diameter), matching the exact r/a grid extracted
  # from beh00-sm01.dat. sort_by='id' preserves the listed order (matches .dat row order).
  [beh00_sm01]
    type = PointValueSampler
    points = '0.0357378 0.0000000 -0.0762000
              0.0333375 0.0000000 -0.0762000
              0.0303657 0.0000000 -0.0762000
              0.0269748 0.0000000 -0.0762000
              0.0238125 0.0000000 -0.0762000
              0.0190500 0.0000000 -0.0762000
              0.0142875 0.0000000 -0.0762000
              0.0095250 0.0000000 -0.0762000
              0.0047625 0.0000000 -0.0762000
              -0.0000000 0.0000000 -0.0762000
              -0.0047625 0.0000000 -0.0762000
              -0.0095250 0.0000000 -0.0762000
              -0.0142875 0.0000000 -0.0762000
              -0.0190500 0.0000000 -0.0762000
              -0.0238125 0.0000000 -0.0762000
              -0.0269748 0.0000000 -0.0762000
              -0.0303657 0.0000000 -0.0762000
              -0.0333375 0.0000000 -0.0762000
              -0.0357378 0.0000000 -0.0762000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_th0225]
    type = PointValueSampler
    points = '-0.0046850 0.0000000 0.2032194
              -0.0069026 0.0000000 0.2023008
              -0.0096482 0.0000000 0.2011636
              -0.0127810 0.0000000 0.1998659
              -0.0157026 0.0000000 0.1986558
              -0.0201026 0.0000000 0.1968332
              -0.0245025 0.0000000 0.1950107
              -0.0289025 0.0000000 0.1931882
              -0.0333025 0.0000000 0.1913656
              -0.0377025 0.0000000 0.1895431
              -0.0421024 0.0000000 0.1877206
              -0.0465024 0.0000000 0.1858980
              -0.0509024 0.0000000 0.1840755
              -0.0553024 0.0000000 0.1822530
              -0.0597023 0.0000000 0.1804305
              -0.0626239 0.0000000 0.1792203
              -0.0657567 0.0000000 0.1779227
              -0.0685023 0.0000000 0.1767854
              -0.0707199 0.0000000 0.1758668'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_th0675]
    type = PointValueSampler
    points = '-0.2920806 0.0000000 0.4906150
              -0.2929992 0.0000000 0.4883974
              -0.2941364 0.0000000 0.4856518
              -0.2954341 0.0000000 0.4825190
              -0.2966442 0.0000000 0.4795974
              -0.2984668 0.0000000 0.4751974
              -0.3002893 0.0000000 0.4707975
              -0.3021118 0.0000000 0.4663975
              -0.3039344 0.0000000 0.4619975
              -0.3057569 0.0000000 0.4575975
              -0.3075794 0.0000000 0.4531976
              -0.3094020 0.0000000 0.4487976
              -0.3112245 0.0000000 0.4443976
              -0.3130470 0.0000000 0.4399976
              -0.3148695 0.0000000 0.4355977
              -0.3160797 0.0000000 0.4326761
              -0.3173773 0.0000000 0.4295433
              -0.3185146 0.0000000 0.4267977
              -0.3194332 0.0000000 0.4245801'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_th1125]
    type = PointValueSampler
    points = '-0.6985194 0.0000000 0.4906150
              -0.6976008 0.0000000 0.4883974
              -0.6964636 0.0000000 0.4856518
              -0.6951659 0.0000000 0.4825190
              -0.6939558 0.0000000 0.4795974
              -0.6921332 0.0000000 0.4751974
              -0.6903107 0.0000000 0.4707975
              -0.6884882 0.0000000 0.4663975
              -0.6866656 0.0000000 0.4619975
              -0.6848431 0.0000000 0.4575975
              -0.6830206 0.0000000 0.4531976
              -0.6811980 0.0000000 0.4487976
              -0.6793755 0.0000000 0.4443976
              -0.6775530 0.0000000 0.4399976
              -0.6757305 0.0000000 0.4355977
              -0.6745203 0.0000000 0.4326761
              -0.6732227 0.0000000 0.4295433
              -0.6720854 0.0000000 0.4267977
              -0.6711668 0.0000000 0.4245801'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_th1575]
    type = PointValueSampler
    points = '-0.9859150 0.0000000 0.2032194
              -0.9836974 0.0000000 0.2023008
              -0.9809518 0.0000000 0.2011636
              -0.9778190 0.0000000 0.1998659
              -0.9748974 0.0000000 0.1986558
              -0.9704974 0.0000000 0.1968332
              -0.9660975 0.0000000 0.1950107
              -0.9616975 0.0000000 0.1931882
              -0.9572975 0.0000000 0.1913656
              -0.9528975 0.0000000 0.1895431
              -0.9484976 0.0000000 0.1877206
              -0.9440976 0.0000000 0.1858980
              -0.9396976 0.0000000 0.1840755
              -0.9352976 0.0000000 0.1822530
              -0.9308977 0.0000000 0.1804305
              -0.9279761 0.0000000 0.1792203
              -0.9248433 0.0000000 0.1779227
              -0.9220977 0.0000000 0.1767854
              -0.9198801 0.0000000 0.1758668'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_sp01]
    type = PointValueSampler
    points = '-1.0263378 0.0000000 -0.0762000
              -1.0239375 0.0000000 -0.0762000
              -1.0209657 0.0000000 -0.0762000
              -1.0175748 0.0000000 -0.0762000
              -1.0144125 0.0000000 -0.0762000
              -1.0096500 0.0000000 -0.0762000
              -1.0048875 0.0000000 -0.0762000
              -1.0001250 0.0000000 -0.0762000
              -0.9953625 0.0000000 -0.0762000
              -0.9906000 0.0000000 -0.0762000
              -0.9858375 0.0000000 -0.0762000
              -0.9810750 0.0000000 -0.0762000
              -0.9763125 0.0000000 -0.0762000
              -0.9715500 0.0000000 -0.0762000
              -0.9667875 0.0000000 -0.0762000
              -0.9636252 0.0000000 -0.0762000
              -0.9602343 0.0000000 -0.0762000
              -0.9572625 0.0000000 -0.0762000
              -0.9548622 0.0000000 -0.0762000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_sp06]
    type = PointValueSampler
    points = '-1.0263378 0.0000000 -0.4572000
              -1.0239375 0.0000000 -0.4572000
              -1.0209657 0.0000000 -0.4572000
              -1.0175748 0.0000000 -0.4572000
              -1.0144125 0.0000000 -0.4572000
              -1.0096500 0.0000000 -0.4572000
              -1.0048875 0.0000000 -0.4572000
              -1.0001250 0.0000000 -0.4572000
              -0.9953625 0.0000000 -0.4572000
              -0.9906000 0.0000000 -0.4572000
              -0.9858375 0.0000000 -0.4572000
              -0.9810750 0.0000000 -0.4572000
              -0.9763125 0.0000000 -0.4572000
              -0.9715500 0.0000000 -0.4572000
              -0.9667875 0.0000000 -0.4572000
              -0.9636252 0.0000000 -0.4572000
              -0.9602343 0.0000000 -0.4572000
              -0.9572625 0.0000000 -0.4572000
              -0.9548622 0.0000000 -0.4572000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_sp10]
    type = PointValueSampler
    points = '-1.0263378 0.0000000 -0.7620000
              -1.0239375 0.0000000 -0.7620000
              -1.0209657 0.0000000 -0.7620000
              -1.0175748 0.0000000 -0.7620000
              -1.0144125 0.0000000 -0.7620000
              -1.0096500 0.0000000 -0.7620000
              -1.0048875 0.0000000 -0.7620000
              -1.0001250 0.0000000 -0.7620000
              -0.9953625 0.0000000 -0.7620000
              -0.9906000 0.0000000 -0.7620000
              -0.9858375 0.0000000 -0.7620000
              -0.9810750 0.0000000 -0.7620000
              -0.9763125 0.0000000 -0.7620000
              -0.9715500 0.0000000 -0.7620000
              -0.9667875 0.0000000 -0.7620000
              -0.9636252 0.0000000 -0.7620000
              -0.9602343 0.0000000 -0.7620000
              -0.9572625 0.0000000 -0.7620000
              -0.9548622 0.0000000 -0.7620000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_sp18]
    type = PointValueSampler
    points = '-1.0263378 0.0000000 -1.3716000
              -1.0239375 0.0000000 -1.3716000
              -1.0209657 0.0000000 -1.3716000
              -1.0175748 0.0000000 -1.3716000
              -1.0144125 0.0000000 -1.3716000
              -1.0096500 0.0000000 -1.3716000
              -1.0048875 0.0000000 -1.3716000
              -1.0001250 0.0000000 -1.3716000
              -0.9953625 0.0000000 -1.3716000
              -0.9906000 0.0000000 -1.3716000
              -0.9858375 0.0000000 -1.3716000
              -0.9810750 0.0000000 -1.3716000
              -0.9763125 0.0000000 -1.3716000
              -0.9715500 0.0000000 -1.3716000
              -0.9667875 0.0000000 -1.3716000
              -0.9636252 0.0000000 -1.3716000
              -0.9602343 0.0000000 -1.3716000
              -0.9572625 0.0000000 -1.3716000
              -0.9548622 0.0000000 -1.3716000'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []


  [beh00_th0675_wall]
    type = PointValueSampler
    points = '-0.3201549 0.0000000 0.4228377
              -0.2913589 0.0000000 0.4923573'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []

  [beh00_th1125_wall]
    type = PointValueSampler
    points = '-0.6704451 0.0000000 0.4228377
              -0.6992411 0.0000000 0.4923573'
    variable = 'vel_x vel_y vel_z pressure TKE TKED yplus'
    sort_by = 'id'
    execute_on = 'FINAL'
  []
  
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
