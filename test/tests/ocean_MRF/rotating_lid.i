### Thermophysical Properties ###
mu = 1.0
rho = 1.0

### Operation Conditions ###
# lid_velocity = 1.0
side_length = 1.0

#walls = 'left top right bottom'

w = 1.0


[GlobalParams]
  #rhie_chow_user_object = 'rc'
  advected_interp_method = 'upwind'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = ${fparse -side_length}
    xmax = ${fparse side_length}
    ymin = ${fparse -side_length}
    ymax = ${fparse side_length}
    nx = 80
    ny = 80
  []
  # Prevent test diffing on distributed parallel element numbering
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
    body_force_kernel_names = "u_omega; v_omega"
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
    initial_condition = 0
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
    momentum_component = 'x'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
  []
  [u_omega]
    type = LinearFVSource
    variable = vel_x
    source_density = x_accel
  []

  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    mu = ${mu}
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
  []
  [v_omega]
    type = LinearFVSource
    variable = vel_y
    source_density = y_accel
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
    boundary = 'left right top bottom'
    functor = 0
  []
  # [pressure-extrapolation]
  #   type = LinearFVExtrapolatedPressureBC
  #   boundary = 'left right top bottom'
  #   variable = pressure
  #   use_two_term_expansion = true
  # []
  [pressure]
    type = LinearFVPressureFluxBC
    boundary = 'top bottom left right'
    variable = pressure
    HbyA_flux = HbyA
    Ainv = Ainv
  []
[]

[AuxVariables]
[]

[AuxKernels]
[]

[FunctorMaterials]
  [x_acceleration]
    type = ParsedFunctorMaterial
    expression = '${w}*${w}*x'
    property_name = 'x_accel'
  []
  [y_acceleration]
    type = ParsedFunctorMaterial
    expression = '${w}*${w}*y'
    property_name = 'y_accel'
  []
[]

[Executioner]
  type = SIMPLE

  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'

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
[]

[Outputs]
  csv = true
  exodus = true
  perf_graph = false
  print_nonlinear_residuals = false
  print_linear_residuals = true
[]

[VectorPostprocessors]
[]
