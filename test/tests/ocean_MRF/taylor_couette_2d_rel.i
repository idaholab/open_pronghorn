### Thermophysical Properties ###
mu = 1.0e-4
rho = 1.0
omega_vel = 0.010016 #factor correcting vel bc at the top due to centroid difference with wall
omega = -0.01

#walls = 'left top right bottom'

[GlobalParams]
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  advected_interp_method = 'upwind'
  u = vel_x
  v = vel_y
[]

[Mesh]
  # [./ccmg]
  #   type = ConcentricCircleMeshGenerator
  #   num_sectors = 400
  #   radii = '0.35 1'
  #   rings = '1 120'
  #   has_outer_square = off
  #   preserve_volumes = off
  #   smoothing_max_it = 3
  # []
  # [add_internal]
  #   type = SideSetsBetweenSubdomainsGenerator
  #   input = ccmg
  #   paired_block = 1
  #   primary_block = 2
  #   new_boundary = inner
  # []
  # [ed0]
  #   type = BlockDeletionGenerator
  #   input = add_internal
  #   block = '1'
  # []
  [file_mesh]
    type = FileMeshGenerator
    file = taylor_couette_2d_rel_out_visc_fine.e
    use_for_exodus_restart = true
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system'
  previous_nl_solution_required = true
[]

[UserObjects]
  [ins_rhie_chow_interpolator]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    # body_force_kernel_names = "u_srfacccel; v_srfacccel"
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    # initial_condition = 0.0
    initial_from_file_var = vel_x
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    # initial_condition = 0
    initial_from_file_var = vel_y
    solver_sys = v_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    # initial_condition = 1e-8
    solver_sys = pressure_system
    initial_from_file_var = pressure
  []
[]

[AuxVariables]
  [wr_0]
    type = MooseLinearVariableFVReal
  []
  [wr_1]
    type = MooseLinearVariableFVReal
  []
  [vel_abs_x]
    type = MooseLinearVariableFVReal
  []
  [vel_abs_y]
    type = MooseLinearVariableFVReal
  []
[]

[AuxKernels]
  [wr_0]
    type = ParsedAux
    variable = wr_0
    use_xyzt = true
    expression = '-y * ${omega_vel}' #wxr=(-wy,wx,0)
    execute_on = 'INITIAL NONLINEAR'
  []
  [wr_1]
    type = ParsedAux
    variable = wr_1
    use_xyzt = true
    expression = 'x * ${omega_vel}'
    execute_on = 'INITIAL NONLINEAR'
  []
  [vel_abs_x]
    type = ParsedAux
    variable = vel_abs_x
    use_xyzt = true
    coupled_variables = 'vel_x'
    expression = 'vel_x - y * ${omega}' #v_abs=v_rel+wxr -->
    execute_on = 'INITIAL NONLINEAR'
  []
  [vel_abs_y]
    type = ParsedAux
    variable = vel_abs_y
    use_xyzt = true
    coupled_variables = 'vel_y'
    expression = 'vel_y + x * ${omega}'
    execute_on = 'INITIAL NONLINEAR'
  []
[]

[LinearFVKernels]
  # [u_time]
  #   type = LinearFVTimeDerivative
  #   variable = vel_x
  #   factor = ${fparse rho}
  # []
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
    use_nonorthogonal_correction = false #true
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
  []
  [u_srfacccel]
    type = LinearFVSRFAccelerations
    variable = vel_x
    omega_brf = omega_brf
    omega_dot_brf = omega_dot_brf
    r_mc = r_mc
    momentum_component = 'x'
    rho = ${rho}
    u = vel_x
    v = vel_y
  []

  # [v_time]
  #   type = LinearFVTimeDerivative
  #   variable = vel_y
  #   factor = ${fparse rho}
  # []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
    use_nonorthogonal_correction = false #true
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
  []
  [v_srfacccel]
    type = LinearFVSRFAccelerations
    variable = vel_y
    omega_brf = omega_brf
    omega_dot_brf = omega_dot_brf
    r_mc = r_mc
    momentum_component = 'y'
    rho = ${rho}
    u = vel_x
    v = vel_y
  []

  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false #true
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = false
  []
[]

[LinearFVBCs]
  [inner-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_x
    boundary = 'inner'
    functor = 0
  []
  [inner-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_y
    boundary = 'inner'
    functor = 0
  []
  [outer-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_x
    boundary = 'outer'
    functor = 'wr_0'
  []
  [outer-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_y
    boundary = 'outer'
    functor = 'wr_1'
  []

  # [pressure]
  #   type = LinearFVAdvectionDiffusionExtrapolatedBC
  #   variable = pressure
  #   boundary = 'inner outer'
  #   use_two_term_expansion = false
  # []
  [pressure]
    type = LinearFVPressureFluxBC
    boundary = 'inner outer'
    variable = pressure
    u = vel_x
    v = vel_y
    rho = ${rho}
    HbyA_flux = HbyA
    Ainv = Ainv
  []
[]

[FunctorMaterials]
  [SRF_Functor_Material]
    type = LinearFVSRFFunctorMaterial
    mc_origin = '0. 0. 0.'
    SRF_input_mode = 'fixed'
    pitch_angle_fixed = 0.0
    roll_angle_fixed = 0.0
    yaw_angle_fixed = 0.0
    pitch_omega_fixed = 0.0
    roll_omega_fixed = 0.0
    yaw_omega_fixed = ${omega}
  []
[]

################################################################################
# EXECUTION / SOLVE
################################################################################

[Executioner]
  # type = PIMPLE
  # momentum_l_abs_tol = 1e-11
  # pressure_l_abs_tol = 1e-11
  # momentum_l_tol = 0
  # pressure_l_tol = 0
  # rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  # momentum_systems = 'u_system v_system'
  # pressure_system = 'pressure_system'
  # momentum_equation_relaxation = 0.7
  # pressure_variable_relaxation = 0.3
  # num_iterations = 6
  # dt = 15
  # num_steps = 2000
  # pressure_absolute_tolerance = 1e-8
  # momentum_absolute_tolerance = 1e-8
  # print_fields = false
  # momentum_l_max_its = 300

  type = SIMPLE
  momentum_l_abs_tol = 1e-18 #1e-11
  pressure_l_abs_tol = 1e-18#1e-11
  momentum_l_tol = 1e-18#0
  pressure_l_tol = 1e-18#0
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.25
  num_iterations = 10000 #6
  #dt = 15
  #num_steps = 2000
  pressure_absolute_tolerance = 1e-16
  momentum_absolute_tolerance = 1e-16
  print_fields = false
  momentum_l_max_its = 300

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.75 0.0 0.0'

  # momentum_petsc_options = '-ksp_monitor'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'

  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'

  continue_on_max_its = true
[]

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  [out_visc]
    type = Exodus
  []
[]
