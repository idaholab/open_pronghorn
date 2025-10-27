### Avg. Fluid Properties for Air at 15-35 degrees Celsius
mu = 1.8485e-5 # Dynamic viscosity [kg/(m*s)]
rho = 1.185 # Density
k = 0.02605 # Thermal conductivity [W/(m*K)]
cp = 1005 # Specific heat cp [J/(kg*K)]
beta = 0.00336 # Fluid expansion coefficient

### Operating Conditions
cold_temp = 288.25 # [Kelvin] Given: 15.1 Celsius
hot_temp = 307.85 # [Kelvin] Given: 34.7 Celsius
pressure_air = 101325 # [Pascals]
length_cavity = 0.076 # [Meters]; Also the reference length used in k-epsilon
u_ref = 0.001 # [m/s]; A reference flow speed
Pr_t = 0.9
# Rayleigh number = 8.6e5

### Initial Conditions
intensity = 0.01 # Turbulence intensity
k_init = '${fparse 1.5 * (intensity * u_ref) ^ 2}' # Turbulent initial kinetic energy
eps_init = '${fparse ((C_mu ^ 0.75) * (k_init ^ 1.5)) / length_cavity}' # Turbulent initial dissipation rate

### k-epsilon Closure Parameters: Widely accepted values from Launder and Sharma (1974)
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09

### Wall Conditions
walls = 'right left top bottom'
bulk_wall_treatment = false
wall_treatment = 'eq_newton'

advected_interp_method = 'upwind'

[Mesh]
  [mesh_cavity]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${length_cavity}
    ymin = 0
    ymax = 2.18
    nx = 10
    ny = 25
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system energy_system TKE_system TKED_system'
  previous_nl_solution_required = true
[]

[Physics]
  [NavierStokes]
    [FlowSegregated]
      [flow]
        compressibility = 'weakly-compressible'
        verbose = true

        velocity_variable = 'velocity_x velocity_y'

        density = 'rho'
        dynamic_viscosity = 'mu'

        initial_velocity = '1e-15 1e-15 0'
        initial_pressure = '${pressure_air}'

        wall_boundaries = 'right left top bottom'
        momentum_wall_types = 'noslip noslip noslip noslip'

        # mass_advection_interpolation = ${advected_interp_method}
        momentum_advection_interpolation = ${advected_interp_method}

        gravity = '0 -9.81 0'
      []
    []
    [FluidHeatTransferSegregated]
      [energy]
        coupled_flow_physics = flow
        verbose = true

        thermal_conductivity = '${k}'
        specific_heat = '${cp}'

        initial_temperature = ${cold_temp}

        energy_wall_types = 'fixed-temperature fixed-temperature heatflux heatflux'
        energy_wall_functors = '${hot_temp} ${cold_temp} 0 0'

        energy_advection_interpolation = ${advected_interp_method}
      []
    []
  []
[]

[Variables]
  [TKE]
  type = MooseLinearVariableFVReal
  solver_sys = TKE_system
  initial_condition = ${k_init}
  []
  [TKED]
  type = MooseLinearVariableFVReal
  solver_sys = TKED_system
  initial_condition = ${eps_init}
  []
[]

[LinearFVKernels]
  [TKE_advection]
    type = LinearFVTurbulentAdvection
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
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
    u = velocity_x
    v = velocity_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 10

    # The parameters below initialize the C3-buoyancy term
    temperature = ${cold_temp}
    alpha_name = ${beta}
    Pr_t = ${Pr_t}
    gravity = '0 -9.81 0'
  []

  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
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
    u = velocity_x
    v = velocity_y
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 10

    # The parameters below initialize the C3-buoyancy term
    temperature = ${cold_temp}
    alpha_name = ${beta}
    Pr_t = ${Pr_t}
    gravity = '0 -9.81 0'
  []
[]

# Kernels for momentum and energy
[LinearFVKernels]
  [turbulent_heat_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = T_fluid
    diffusion_coeff = 'k_t'
    use_nonorthogonal_correction = true
    walls = ${walls}
  []
[]

[LinearFVBCs]
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = ${walls}
    variable = 'mu_t'
    u = velocity_x
    v = velocity_y
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
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
  [k_t]
    type = MooseLinearVariableFVReal
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
    u = velocity_x
    v = velocity_y
    bulk_wall_treatment = ${bulk_wall_treatment}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
    mu_t_ratio_max = 1e20
  []
  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    tke = TKE
    mu = ${mu}
    rho = ${rho}
    u = velocity_x
    v = velocity_y
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
  [compute_k_t]
    type = TurbulentConductivityAux
    cp = ${cp}
    mu_t = 'mu_t'
    variable = 'k_t'
    Pr_t = ${Pr_t}
  []
[]

[FunctorMaterials]
  [const_functor]
    type = ADGenericFunctorMaterial
    prop_names = 'cp    k  mu  drho_dt'
    prop_values = '${cp} ${k} ${mu} 0'
  []
  [density]
    type = ADParsedFunctorMaterial
    property_name = 'rho'
    expression = '${rho} * (1 - ${beta} * (T_fluid - ${cold_temp}))'
    functor_names = 'T_fluid'
  []
  [density_cp]
    type = ADParsedFunctorMaterial
    property_name = 'rho_cp'
    expression = '${cp} * ${rho} * (1 - ${beta} * (T_fluid - ${cold_temp}))'
    functor_names = 'T_fluid'
  []
[]

[Executioner]
  type = PIMPLE

  momentum_l_abs_tol = 1e-11
  pressure_l_abs_tol = 1e-11
  energy_l_abs_tol = 1e-11
  momentum_l_tol = 1e-8
  pressure_l_tol = 1e-8
  energy_l_tol = 1e-8
  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  energy_absolute_tolerance = 1e-8

  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  energy_system = 'energy_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  energy_equation_relaxation = 0.9

  num_iterations = 10
  print_fields = false
  continue_on_max_its = true
  dt = 0.01
  num_steps = 5000
  num_piso_iterations = 0
#   momentum_l_max_its = 300

  pin_pressure = true
  pressure_pin_value = 101325
  pressure_pin_point = '0.038 1.09 0.0'

  # momentum_petsc_options = '-ksp_monitor'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'

  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'

  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
[]

[VectorPostprocessors]

  ### .csv Files for y/H = 0.10
  [plot_temperature_10]
    type = LineValueSampler
    variable = T_fluid
    start_point = '0 0.218 0'
    end_point = '0.076 0.218 0'
    num_points = 25
    sort_by = x
  []
  [plot_vel_x_10]
    type = LineValueSampler
    variable = velocity_x
    start_point = '0 0.218 0'
    end_point = '0.076 0.218 0'
    num_points = 25
    sort_by = x
  []
  [plot_vel_y_10]
    type = LineValueSampler
    variable = velocity_y
    start_point = '0 0.218 0'
    end_point = '0.076 0.218 0'
    num_points = 25
    sort_by = x
  []

  ### .csv Files for y/H = 0.50
  [plot_temperature_50]
    type = LineValueSampler
    variable = T_fluid
    start_point = '0 1.09 0'
    end_point = '0.076 1.09 0'
    num_points = 25
    sort_by = x
  []
  [plot_vel_x_50]
    type = LineValueSampler
    variable = velocity_x
    start_point = '0 1.09 0'
    end_point = '0.076 1.09 0'
    num_points = 25
    sort_by = x
  []
  [plot_vel_y_50]
    type = LineValueSampler
    variable = velocity_y
    start_point = '0 1.09 0'
    end_point = '0.076 1.09 0'
    num_points = 25
    sort_by = x
  []

  ### .csv Files for y/H = 0.90
  [plot_temperature_90]
    type = LineValueSampler
    variable = T_fluid
    start_point = '0 1.962 0'
    end_point = '0.076 1.962 0'
    num_points = 25
    sort_by = x
  []
  [plot_vel_x_90]
    type = LineValueSampler
    variable = velocity_x
    start_point = '0 1.962 0'
    end_point = '0.076 1.962 0'
    num_points = 25
    sort_by = x
  []
  [plot_vel_y_90]
    type = LineValueSampler
    variable = velocity_y
    start_point = '0 1.962 0'
    end_point = '0.076 1.962 0'
    num_points = 25
    sort_by = x
  []
[]


[Outputs]
  exodus = true

  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
