################################################################################
# MATERIAL PROPERTIES
################################################################################
rho_0 = 3279.
T_0 = 875.0
mu = 0.001
k_cond = 38.0
cp = 640.
beta = 3.26e-5

walls = 'right left top bottom'

[GlobalParams]
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  advected_interp_method = 'upwind'
  u = vel_x
  v = vel_y
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system energy_system'
  previous_nl_solution_required = true
[]

################################################################################
# GEOMETRY
################################################################################

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = 100
    ny = 100
  []
[]

################################################################################
# EQUATIONS: VARIABLES, KERNELS & BCS
################################################################################

[UserObjects]
  [ins_rhie_chow_interpolator]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho_0}
    p_diffusion_kernel = p_diffusion
    # body_force_kernel_names = "u_buoyancy u_srfacccel; v_buoyancy v_srfacccel"
    #body_force_kernel_names = "u_srfacccel; v_srfacccel"
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
    initial_condition = 0
    solver_sys = pressure_system
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_condition = 875
  []
[]

[LinearFVKernels]
  [u_time]
    type = LinearFVTimeDerivative
    variable = vel_x
    factor = ${fparse rho_0}
  []
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
    use_nonorthogonal_correction = false
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
  []
  [u_buoyancy]
    type = LinearFVSRFMomentumBoussinesq
    variable = vel_x
    T_fluid = T_fluid
    gravity = '0 0.0 0'
    rho = ${rho_0}
    ref_temperature = ${T_0}
    alpha_name = ${beta}
    momentum_component = 'x'
    pitch_angle = pitch_angle
    roll_angle = roll_angle
    yaw_angle = yaw_angle
  []
  [u_srfacccel]
    type = LinearFVSRFAccelerations
    variable = vel_x
    omega_brf = omega_brf
    omega_dot_brf = omega_dot_brf
    r_mc = r_mc
    momentum_component = 'x'
    rho = ${rho_0}
    u = vel_x
    v = vel_y
  []

  [v_time]
    type = LinearFVTimeDerivative
    variable = vel_y
    factor = ${fparse rho_0}
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
    use_nonorthogonal_correction = false
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
  []
  [v_buoyancy]
    type = LinearFVSRFMomentumBoussinesq
    variable = vel_y
    T_fluid = T_fluid
    gravity = '0 0.0 0'
    rho = ${rho_0}
    ref_temperature = ${T_0}
    alpha_name = ${beta}
    momentum_component = 'y'
    pitch_angle = pitch_angle
    roll_angle = roll_angle
    yaw_angle = yaw_angle
  []
  [v_srfacccel]
    type = LinearFVSRFAccelerations
    variable = vel_y
    omega_brf = omega_brf
    omega_dot_brf = omega_dot_brf
    r_mc = r_mc
    momentum_component = 'y'
    rho = ${rho_0}
    u = vel_x
    v = vel_y
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
    force_boundary_execution = false
  []

   ####### FUEL ENERGY EQUATION #######
  [T_time]
    type = LinearFVTimeDerivative
    variable = T_fluid
    factor = ${fparse rho_0*cp}
  []
  [heat_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
  []
  [conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = ${fparse k_cond}
  []
[]


[LinearFVBCs]
  [no-slip-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_x
    boundary = ${walls}
    functor = 0
  []
  [no-slip-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_y
    boundary = ${walls}
    functor = 0
  []
  [T_cold]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = T_fluid
    boundary = 'right'
    functor = 870.0
  []
  [T_hot]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = T_fluid
    boundary = 'left'
    functor = 880.0
  []
  [T_all]
    type = LinearFVAdvectionDiffusionExtrapolatedBC
    variable = T_fluid
    boundary = 'top bottom'
    use_two_term_expansion = false
  []
  # [pressure]
  #   type = LinearFVAdvectionDiffusionExtrapolatedBC
  #   variable = pressure
  #   boundary = 'top bottom left right'
  #   use_two_term_expansion = false
  # []
  [pressure]
    type = LinearFVPressureFluxBC
    boundary = 'top bottom left right'
    variable = pressure
    u = vel_x
    v = vel_y
    rho = ${rho_0}
    HbyA_flux = HbyA
    Ainv = Ainv
  []
[]

[FunctorMaterials]
  [constant_functors]
    type = GenericFunctorMaterial
    prop_names = 'cp beta'
    prop_values = '${cp} ${beta}'
  []
  [SRF_Functor_Material]
    type = LinearFVSRFFunctorMaterial
    mc_origin = '0.5 0.5 0'
    SRF_input_mode = 'fixed'
    pitch_angle_fixed = 0.0
    roll_angle_fixed = 0.0
    yaw_angle_fixed = 0.0
    pitch_omega_fixed = 0.0
    roll_omega_fixed = 0.0
    yaw_omega_fixed = 0.1
  []
[]

################################################################################
# EXECUTION / SOLVE
################################################################################

[Executioner]
  type = PIMPLE
  momentum_l_abs_tol = 1e-11
  pressure_l_abs_tol = 1e-11
  energy_l_abs_tol = 1e-11
  momentum_l_tol = 0
  pressure_l_tol = 0
  energy_l_tol = 0
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  energy_system = 'energy_system'
  momentum_equation_relaxation = 0.6
  pressure_variable_relaxation = 0.25
  energy_equation_relaxation = 0.8
  num_iterations = 10 #8000
  dt = 1
  num_steps = 3000
  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  energy_absolute_tolerance = 1e-8
  print_fields = false
  momentum_l_max_its = 300

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.5 0.5 0.0'

  # momentum_petsc_options = '-ksp_monitor'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'

  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'

  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'

  continue_on_max_its = true
[]

################################################################################
# SIMULATION OUTPUTS
################################################################################

[Outputs]
  //exodus = true
  [rot]
    type = Exodus
  []
[]
