# Fixed-angular-velocity SRF equilibrium.
#
# For omega = (0, 0, w), mc_origin = (0, 0, 0), and zero relative
# velocity, the centrifugal acceleration is
#
#   -omega x (omega x r) = w^2 (x, y, 0).
#
# The steady analytical solution used by the error postprocessors is
#
#   u = 0,
#   v = 0,
#   p = 0.5 * rho * w^2 * (x^2 + y^2),
#
# where the pressure constant is selected by pinning p(0, 0) = 0.

mu = 1.0
rho = 1.0
side_length = 1.0
w = 1.0

[FVInterpolationMethods]
  [average]
    type = FVGeometricAverage
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${fparse -side_length}
    xmax = ${fparse side_length}
    ymin = ${fparse -side_length}
    ymax = ${fparse side_length}
    # mms.run_spatial uniformly refines this 2 x 2 base mesh.
    nx = 2
    ny = 2
  []
  # Prevent test diffing on distributed parallel element numbering.
  allow_renumbering = false
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system'
  previous_nl_solution_required = true
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
    initial_condition = 0.0
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    initial_condition = 0.0
    solver_sys = v_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    initial_condition = 1e-8
    solver_sys = pressure_system
  []
[]

[LinearFVKernels]
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = ${mu}
    u = vel_x
    v = vel_y
    momentum_component = x
    rhie_chow_user_object = rc
    use_nonorthogonal_correction = false
    use_deviatoric_terms = true
    advected_interp_method_name = average
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = x
  []
  [u_omega]
    type = LinearFVSRFAccelerations
    variable = vel_x
    momentum_component = x
    rho = ${rho}
    u = vel_x
    v = vel_y
    omega_brf = omega_brf
    omega_dot_brf = omega_dot_brf
    r_mc = r_mc
  []

  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = ${mu}
    u = vel_x
    v = vel_y
    momentum_component = y
    rhie_chow_user_object = rc
    use_nonorthogonal_correction = false
    use_deviatoric_terms = true
    advected_interp_method_name = average
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = y
  []
  [v_omega]
    type = LinearFVSRFAccelerations
    variable = vel_y
    momentum_component = y
    rho = ${rho}
    u = vel_x
    v = vel_y
    omega_brf = omega_brf
    omega_dot_brf = omega_dot_brf
    r_mc = r_mc
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
  [no_slip_x]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_x
    boundary = 'left right bottom top'
    functor = 0
  []
  [no_slip_y]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_y
    boundary = 'left right bottom top'
    functor = 0
  []
  [pressure]
    type = LinearFVPressureFluxBC
    boundary = 'left right bottom top'
    variable = pressure
    HbyA_flux = HbyA
    Ainv = Ainv
    u = vel_x
    v = vel_y
    rho = ${rho}
  []
[]

[Functions]
  [exact_u]
    type = ParsedFunction
    expression = '0'
  []
  [exact_v]
    type = ParsedFunction
    expression = '0'
  []
  [exact_p]
    type = ParsedFunction
    expression = '0.5*rho*w^2*(x^2+y^2)'
    symbol_names = 'rho w'
    symbol_values = '${rho} ${w}'
  []
[]

[FunctorMaterials]
  [srf_motion]
    type = LinearFVSRFFunctorMaterial
    mc_origin = '0 0 0'
    SRF_input_mode = fixed

    pitch_angle_fixed = 0
    yaw_angle_fixed = 0
    roll_angle_fixed = 0

    pitch_omega_fixed = 0
    yaw_omega_fixed = ${w}
    roll_omega_fixed = 0

    pitch_omegadot_fixed = 0
    yaw_omegadot_fixed = 0
    roll_omegadot_fixed = 0
  []
[]

[Executioner]
  type = SIMPLE

  rhie_chow_user_object = rc
  momentum_systems = 'u_system v_system'
  pressure_system = pressure_system

  momentum_l_abs_tol = 1e-14
  pressure_l_abs_tol = 1e-14
  momentum_l_tol = 1e-14
  pressure_l_tol = 1e-14

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  num_iterations = 1000
  pressure_absolute_tolerance = 1e-12
  momentum_absolute_tolerance = 1e-12

  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'

  print_fields = false
  continue_on_max_its = true

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.0 0.0 0.0'

  # Accurately integrate the squared quadratic pressure error.
  [Quadrature]
    type = GAUSS
    order = FOURTH
  []
[]

[Postprocessors]
  [h]
    type = AverageElementSize
    outputs = csv
  []
  [L2u]
    type = ElementL2FunctorError
    approximate = vel_x
    exact = exact_u
    outputs = csv
  []
  [L2v]
    type = ElementL2FunctorError
    approximate = vel_y
    exact = exact_v
    outputs = csv
  []
  [L2p]
    type = ElementL2FunctorError
    approximate = pressure
    exact = exact_p
    outputs = csv
  []
[]

[Outputs]
  csv = true
[]

