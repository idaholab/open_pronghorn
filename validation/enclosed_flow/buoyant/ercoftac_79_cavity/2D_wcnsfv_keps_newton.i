# ERCOFTAC Case No. 079: Natural Convection of Air in a Differentially-Heated Enclosed Cavity

### Avg. Fluid Properties for Air at 15-35 degrees Celsius
mu = 1.8485e-5 # Dynamic viscosity [kg/(m*s)]
rho = 1.185 # Density
k = 0.02605 # Thermal conductivity [W/(m*K)]
cp = 1005 # Specific heat cp [J/(kg*K)]
alpha = 0.00336 # Fluid expansion coefficient

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
    nx = 50
    ny = 50
  []
  uniform_refine = 0
[]

[AuxVariables]
  [yplus]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
  []
  [Reynolds]
    order = CONSTANT # These are the same as type = MooseVariableFVReal
    family = MONOMIAL
    fv = true
  []
[]

[AuxKernels]
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
  [Reynolds]
    type = ReynoldsNumberFunctorAux
    variable = Reynolds
    speed = speed
    rho = ${rho}
    mu = ${mu}
  []
[]

[Physics]
  [NavierStokes]
    [Flow] # FlowSegregated
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

        mass_advection_interpolation = ${advected_interp_method}
        momentum_advection_interpolation = ${advected_interp_method}

        gravity = '0 -9.81 0'
      []
    []
    [FluidHeatTransfer]
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
    [Turbulence]
      [k_epsilon]
        verbose = true

        turbulence_handling = 'k-epsilon'

        transient = true

        tke_name = TKE
        tked_name = TKED

        # Initialization
        initial_tke = ${k_init}
        initial_tked = ${eps_init}

        # Model parameters
        C1_eps = ${C1_eps}
        C2_eps = ${C2_eps}
        C_mu = ${C_mu}

        sigma_k = ${sigma_k}
        sigma_eps = ${sigma_eps}

        # Wall parameters
        turbulence_walls = ${walls}
        bulk_wall_treatment = ${bulk_wall_treatment}
        wall_treatment_eps = ${wall_treatment}
        wall_treatment_T = ${wall_treatment} # ADDED 6/15/2025

        # Numerical parameters
        turbulent_viscosity_two_term_bc_expansion = false
        turbulent_viscosity_interp_method = 'average'
        mu_t_as_aux_variable = false
        output_mu_t = false

        coupled_flow_physics = flow
        fluid_heat_transfer_physics = energy
      []
    []
  []
[]

# [FVBCs]
#   [walls_TKE]
#   type = FVDirichletBC
#   boundary = 'right left top bottom'
#   variable = TKE
#   value = ${k_init}
#   []
#   [walls_TKED]
#   type = FVDirichletBC
#   boundary = 'right left top bottom'
#   variable = TKED
#   value = ${eps_init}
#   []
# []

[FunctorMaterials]
  [const_functor]
    type = ADGenericFunctorMaterial
    prop_names = 'cp    k  mu  drho_dt alpha'
    prop_values = '${cp} ${k} ${mu} 0 ${alpha}'
  []
  [density]
    type = ADParsedFunctorMaterial
    property_name = 'rho'
    expression = '${rho} * (1 - ${alpha} * (T_fluid - ${cold_temp}))'
    functor_names = 'T_fluid'
  []
  [Prandtl]
    type = ADGenericFunctorMaterial
    prop_names = 'Pr_t'
    prop_values = ${Pr_t}
  []
[]

[Problem]
  error_on_jacobian_nonzero_reallocation = false # Resolves error: Use MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE) to turn off this check
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu     NONZERO'

  # Try to use a better time stepper
  [TimeStepper]
    type = SolutionTimeAdaptiveDT
    #type = IterationAdaptiveDT
    dt = 0.3375
    #optimal_iterations = 6
  []
  steady_state_detection = true

  nl_abs_tol = 1e-7
  nl_max_its = 25
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  # compute_scaling_once = false
[]

[Postprocessors]
  [rayleigh_number]
    type = RayleighNumber
    T_cold = ${cold_temp}
    T_hot = ${hot_temp}
    rho_ave = ${rho}
    beta = ${alpha}
    l = ${length_cavity}
    mu_ave = ${mu}
    k_ave = ${k}
    cp_ave = ${cp}
    gravity_magnitude = 9.81
  []
  [T_min]
    type = ADElementExtremeFunctorValue
    functor = 'T_fluid'
    value_type = 'min'
  []
  [T_max]
    type = ADElementExtremeFunctorValue
    functor = 'T_fluid'
    value_type = 'max'
  []
[]

[VectorPostprocessors]

  ### .csv Files for y/H = 0.10
  [plot_temperature_10]
    type = LineValueSampler
    variable = T_fluid
    start_point = '0 0.218 0'
    end_point = '0.076 0.218 0'
    num_points = 100
    sort_by = x
  []
  [plot_vel_x_10]
    type = LineValueSampler
    variable = velocity_x
    start_point = '0 0.218 0'
    end_point = '0.076 0.218 0'
    num_points = 100
    sort_by = x
  []
  [plot_vel_y_10]
    type = LineValueSampler
    variable = velocity_y
    start_point = '0 0.218 0'
    end_point = '0.076 0.218 0'
    num_points = 100
    sort_by = x
  []

  ### .csv Files for y/H = 0.50
  [plot_temperature_50]
    type = LineValueSampler
    variable = T_fluid
    start_point = '0 1.09 0'
    end_point = '0.076 1.09 0'
    num_points = 100
    sort_by = x
  []
  [plot_vel_x_50]
    type = LineValueSampler
    variable = velocity_x
    start_point = '0 1.09 0'
    end_point = '0.076 1.09 0'
    num_points = 100
    sort_by = x
  []
  [plot_vel_y_50]
    type = LineValueSampler
    variable = velocity_y
    start_point = '0 1.09 0'
    end_point = '0.076 1.09 0'
    num_points = 100
    sort_by = x
  []

  ### .csv Files for y/H = 0.90
  [plot_temperature_90]
    type = LineValueSampler
    variable = T_fluid
    start_point = '0 1.962 0'
    end_point = '0.076 1.962 0'
    num_points = 100
    sort_by = x
  []
  [plot_vel_x_90]
    type = LineValueSampler
    variable = velocity_x
    start_point = '0 1.962 0'
    end_point = '0.076 1.962 0'
    num_points = 100
    sort_by = x
  []
  [plot_vel_y_90]
    type = LineValueSampler
    variable = velocity_y
    start_point = '0 1.962 0'
    end_point = '0.076 1.962 0'
    num_points = 100
    sort_by = x
  []
[]

[Outputs]
  exodus = true

  [CSV]
    type = CSV # start_time = time output operates; end_time = time output stops operating
    execute_on = 'FINAL' # time_step_interval = interval (# of time steps) at which output occurs
  []
[]
