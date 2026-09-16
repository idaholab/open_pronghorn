rho = 1.0
bulk_u = 10.
D = 0.0762
mu = 1.524e-5

# advected_interp_method = 'upwind'
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

    parallel_type = distributed
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
    tke = TKE
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
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    w = vel_z
    walls = ${walls}
    bulk_wall_treatment = ${bulk_wall_treatment}
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
    variable = 'vel_x vel_y vel_z pressure TKE TKED'
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
    variable = 'vel_x vel_y vel_z pressure TKE TKED'
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
