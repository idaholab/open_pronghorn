rho_0 = 1.
mu = 1.

[FVInterpolationMethods]
  [upwind]
    type = FVAdvectedUpwind
  []
[]


[GlobalParams]
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  advected_interp_method = 'average'
  u = vel_x
  v = vel_y
  w = vel_z
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = -0.5
    xmax = 0.5
    nx = 10
  []
  # Prevent test diffing on distributed parallel element numbering
  allow_renumbering = false
[]

[Problem]
  linear_sys_names = 'u_system pressure_system'
  previous_nl_solution_required = true
[]

[UserObjects]
  [ins_rhie_chow_interpolator]
    type = RhieChowMassFlux
    u = vel_x
    pressure = pressure
    rho = ${rho_0}
    p_diffusion_kernel = p_diffusion
    #body_force_kernel_names = "u_omega; v_omega"
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    initial_condition = 0.0
    # initial_from_file_var = vel_x
    solver_sys = u_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    initial_condition = 1e-8
    solver_sys = pressure_system
    #initial_from_file_var = pressure
  []
[]

[AuxVariables]
  [roll_angle_var]
    type = MooseVariableFVReal
  []
  [pitch_angle_var]
    type = MooseVariableFVReal
  []
  [yaw_angle_var]
    type = MooseVariableFVReal
  []
  [omega_roll_var]
    type = MooseVariableFVReal
  []
  [omega_pitch_var]
    type = MooseVariableFVReal
  []
  [omega_yaw_var]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [roll_aux]
    type = FunctorAux
    functor = roll_angle
    variable = roll_angle_var
  []
  [pitch_aux]
    type = FunctorAux
    functor = pitch_angle
    variable = pitch_angle_var
  []
  [yaw_aux]
    type = FunctorAux
    functor = yaw_angle
    variable = yaw_angle_var
  []
  [omega_roll_aux]
    type = FunctorAux
    functor = omega_roll
    variable = omega_roll_var
  []
  [omega_pitch_aux]
    type = FunctorAux
    functor = omega_pitch
    variable = omega_pitch_var
  []
  [omega_yaw_aux]
    type = FunctorAux
    functor = omega_yaw
    variable = omega_yaw_var
  []
[]

[LinearFVKernels]
  # [u_time]
  #   type = LinearFVTimeDerivative
  #   variable = vel_x
  #   factor = ${fparse rho_0}
  # []
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    advected_interp_method_name = 'upwind'
  []
  [v_buoyancy]
    type = LinearFVSRFMomentumBoussinesq
    variable = vel_x
    T_fluid = pressure
    gravity = '0 0 -9.81'
    rho = 1.
    ref_temperature = 1.
    alpha_name = 1.
    momentum_component = 'x'
    pitch_angle = pitch_angle
    roll_angle = roll_angle
    yaw_angle = yaw_angle
  []

  [p_diffusion]
    type = LinearFVPressureCorrectionDiffusion
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
[]

[LinearFVBCs]
[]

[Postprocessors]
  [roll]
    type = PointValue
    point = '0. 0. 0'
    variable = roll_angle_var
    execute_on = 'TIMESTEP_END'
  []
  [pitch]
    type = PointValue
    point = '0. 0. 0'
    variable = pitch_angle_var
    execute_on = 'TIMESTEP_END'
  []
  [yaw]
    type = PointValue
    point = '0. 0. 0'
    variable = yaw_angle_var
    execute_on = 'TIMESTEP_END'
  []
  [omega_roll_pp]
    type = PointValue
    point = '0. 0. 0'
    variable = omega_roll_var
    execute_on = 'TIMESTEP_END'
  []
  [omega_pitch_pp]
    type = PointValue
    point = '0. 0. 0'
    variable = omega_pitch_var
    execute_on = 'TIMESTEP_END'
  []
  [omega_yaw_pp]
    type = PointValue
    point = '0. 0. 0'
    variable = omega_yaw_var
    execute_on = 'TIMESTEP_END'
  []
[]

[FunctorMaterials]
  [SRF_Functor_Material]
    type = LinearFVSRFFunctorMaterial
    mc_origin = '0. 0. 0'
    SRF_input_mode = 'pitch_yaw_roll'
    pitch_amp = 10.0
    pitch_per = 15.0
    pitch_pha = -90
    yaw_amp = 5.0
    yaw_per = 20.0
    yaw_pha = 90
    roll_amp = 20.0
    roll_per = 10.0
    roll_pha = 0
  []
[]

################################################################################
# EXECUTION / SOLVE
################################################################################

[Executioner]
  type = PIMPLE
  # should_solve_momentum = false
  # should_solve_pressure = false
  momentum_l_abs_tol = 1e-11
  pressure_l_abs_tol = 1e-11
  momentum_l_tol = 0
  pressure_l_tol = 0

  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  momentum_systems = 'u_system'
  pressure_system = 'pressure_system'
  momentum_equation_relaxation = 0.65
  pressure_variable_relaxation = 0.25
  num_iterations = 1
  dt = 1.0
  num_steps = 12
  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  print_fields = false
  momentum_l_max_its = 300

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.0 0.0 0.0'

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
  csv = true
  exodus = true
[]
