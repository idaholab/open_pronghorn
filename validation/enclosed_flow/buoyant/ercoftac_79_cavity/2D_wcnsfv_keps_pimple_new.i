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
# Rayleigh number = 8.6e5

# Numerical scheme
advected_interp_method = 'upwind'

### Initial Conditions
intensity = 0.01 # Turbulence intensity
k_init = '${fparse 1.5 * (intensity * u_ref) ^ 2}' # Turbulent initial kinetic energy
eps_init = '${fparse ((C_mu ^ 0.75) * (k_init ^ 1.5)) / length_cavity}' # Turbulent initial dissipation rate

### Wall Conditions
walls = 'right left top bottom'
bulk_wall_treatment = false
wall_treatment = 'eq_newton'

# Turbulence-model knobs (optional)
k_epsilon_variant   = 'StandardTwoLayer'  # Standard | StandardLowRe | StandardTwoLayer | Realizable | RealizableTwoLayer
# two_layer_flavor    = 'Wolfstein' # Wolfstein | NorrisReynolds | Xu (only used for *TwoLayer variants)
use_buoyancy        = true
use_compressibility = false
nonlinear_model     = 'none'
curvature_model     = 'none'
use_yap             = false
use_low_re_Gprime   = false

### k-epsilon Closure Parameters
# Widely accepted values from Launder and Sharma (1974)
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09
Pr_t = 0.9

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

        pressure_projection_method = 'consistent'
        solve_for_dynamic_pressure = true
        gravity = '0 -9.81 0'
        # Used by the turbulent buoyancy production in the turbulence Physics
        thermal_expansion = ${beta}
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
    [AdvancedKEpsilon/k-epsilon]
      turbulence_handling = 'k-epsilon'
      tke_name = TKE
      tked_name = TKED
      system_names = 'TKE_system TKED_system'
      fluid_heat_transfer_physics = energy
      verbose = true

      initial_tke = ${k_init}
      initial_tked = ${eps_init}
      initial_mu_t = '${fparse rho * C_mu * k_init^2 / eps_init}'

      # Model parameters
      sigma_k = ${sigma_k}
      sigma_eps = ${sigma_eps}
      C1_eps = ${C1_eps}
      C2_eps = ${C2_eps}
      # 'C_mu' is left at its 0.09 default: naming it here would make the brace substitution
      # resolve against this block's own parameter rather than the global of the same name
      # The two-layer epsilon equation does not apply the production limiter, so the single
      # value forwarded to both equations reproduces the per-kernel 10 / 2 of the kernel setup
      C_pl = 10

      # wall treatment
      turbulence_walls = ${walls}
      wall_treatment_eps = ${wall_treatment}
      bulk_wall_treatment = ${bulk_wall_treatment}

      # k-epsilon variant options. The two-layer formulation needs a wall distance, which the
      # Physics creates on the turbulence walls
      k_epsilon_variant   = ${k_epsilon_variant}
      use_buoyancy        = ${use_buoyancy}
      use_compressibility = ${use_compressibility}
      nonlinear_model     = ${nonlinear_model}
      curvature_model     = ${curvature_model}
      use_yap             = ${use_yap}
      use_low_re_Gprime   = ${use_low_re_Gprime}

      # The buoyancy production takes the temperature from the energy Physics, and the
      # expansion coefficient and gravity vector from the flow Physics
      Pr_t = ${Pr_t}
    []
  []
[]

# The k and epsilon equations are created by the turbulence Physics

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

# The turbulent viscosity wall function is created by the turbulence Physics
# TODO: add a wall function for T_fluid

[AuxVariables]
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
  # The turbulent viscosity and the wall distance are computed by the turbulence Physics
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
    type = EffectiveTurbulentThermalConductivityAux
    cp = ${cp}
    mu_t = 'mu_t'
    k = ${k}
    variable = 'k_t'
    Pr_t = ${Pr_t}
    # The heat diffusion kernel below carries the molecular conductivity separately, so only
    # the turbulent part cp mu_t / Pr_t belongs here
    turbulent_thermal_conductivity = true
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
    # weakly compressible option
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
  turbulence_l_abs_tol = 1e-8
  momentum_l_tol = 1e-8
  pressure_l_tol = 1e-8
  energy_l_tol = 1e-8
  turbulence_l_tol = 1e-10

  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  energy_absolute_tolerance = 1e-8
  turbulence_absolute_tolerance = '1e-8 1e-8'

  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  energy_system = 'energy_system'
  turbulence_systems = 'TKE_system TKED_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  energy_equation_relaxation = 0.9
  turbulence_equation_relaxation = '0.4 0.4'

  num_iterations = 1000
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

  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'

  momentum_l_max_its = 300
  pressure_l_max_its = 300
  turbulence_l_max_its = 30
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
