##########################################################
# ERCOFTAC test case for turbulent channel flow
# Case Number: 032
# Author: Dr. Mauricio Tano & Hailey Tran Kieu
# Last Update: November, 2023 & 2025
# Turbulent model using:
# k-epsilon
# Equilibrium + Newton wall treatment
# SIMPLE solve
##########################################################

H = 1 #halfwidth of the channel
L = 120

Re = 22250

rho = 1
bulk_u = 1
mu = '${fparse rho * bulk_u * 2 * H / Re}'

advected_interp_method = 'upwind'

### k-epsilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09

### Initial and Boundary Conditions ###
intensity = 0.01
k_init = '${fparse 1.5*(intensity * bulk_u)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (2*H)}'

### Modeling parameters ###
bulk_wall_treatment = false
walls = 'bottom top'
wall_treatment = 'eq_newton' # Options: eq_newton, eq_incremental, eq_linearized, neq

[Mesh]
  [block_1]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${L}'
    dy = '0.84 0.16' # '0.75 0.25'
    ix = '50'
    iy = '35 1'
  []
  [block_2_base]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${L}'
    dy = '0.16 0.84' # '0.25 0.75'
    ix = '50'
    iy = '1 35'
  []
  [block_2]
    type = TransformGenerator
    input = block_2_base
    transform = TRANSLATE
    vector_value = '0 -1 0'
  []
  [smg]
    type = StitchedMeshGenerator
    inputs = 'block_1 block_2'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'bottom top'
    merge_boundaries_with_same_name = true
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system'
  previous_nl_solution_required = true
[]

[Physics]
  [NavierStokes]
    [FlowSegregated/flow]
      velocity_variable = 'vel_x vel_y'
      pressure_variable = 'pressure'

      # Initial conditions
      initial_velocity = '${bulk_u} 0 0'
      initial_pressure = '1e-8'

      # Material properties
      density = ${rho}
      dynamic_viscosity = ${mu}

      # Boundary conditions
      inlet_boundaries = 'left'
      momentum_inlet_types = 'fixed-velocity'
      momentum_inlet_functors = '${bulk_u} 0'

      wall_boundaries = 'top bottom'
      momentum_wall_types = 'noslip noslip'

      outlet_boundaries = 'right'
      momentum_outlet_types = 'fixed-pressure'
      pressure_functors = '0'

      # Numerical parameters
      include_deviatoric_stress = true
      orthogonality_correction = false
      pressure_two_term_bc_expansion = false
      momentum_two_term_bc_expansion = false
      momentum_advection_interpolation = ${advected_interp_method}
      system_names = 'u_system v_system pressure_system'
      # consistency with kernel syntax input, allows to use the same executioner block
      rhie_chow_uo_name = 'rc'
    []
    [TurbulenceSegregated/k-epsilon]
      # Model
      turbulence_handling = 'k-epsilon'
      tke_name = TKE
      tked_name = TKED
      system_names = 'TKE_system TKED_system'

      initial_tke = ${k_init}
      initial_tked = ${eps_init}

      # Model parameters
      mu_t_ratio_max = 1e20
      sigma_k = ${sigma_k}
      sigma_eps = ${sigma_eps}
      C_pl = 1e10
      C1_eps = ${C1_eps}
      C2_eps = ${C2_eps}

      turbulence_walls = ${walls}
      wall_treatment_eps = ${wall_treatment}
      bulk_wall_treatment = ${bulk_wall_treatment}
      use_nonorthogonal_correction = false
    []
  []
[]

[LinearFVBCs]
  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = TKE
    functor = '${k_init}'
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'right'
    variable = TKE
    use_two_term_expansion = false
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = TKED
    functor = '${eps_init}'
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'right'
    variable = TKED
    use_two_term_expansion = false
  []
[]

[AuxVariables]
  [yplus]
    type = MooseVariableFVReal
        two_term_boundary_expansion = false
  []
  [mu_eff_var]
    type = MooseVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
    two_term_boundary_expansion = false
  []
[]

[AuxKernels]
  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    tke = TKE
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []
  [compute_mu_eff_var]
    type = ParsedAux
    variable = 'mu_eff_var'
    coupled_variables = 'mu_t'
    expression = 'mu_t + ${mu}'
  []
[]

[VectorPostprocessors]
  [test_csv]
    type = PointValueSampler
    variable = 'mu_eff_var mu_t pressure TKE TKED vel_x vel_y yplus'
    points = '119.95 0.012 0   119.95 0.036 0   119.95 0.06 0   119.95 0.084 0   119.95 0.108 0
              119.95 0.132 0   119.95 0.156 0   119.95 0.18 0   119.95 0.204 0   119.95 0.228 0
              119.95 0.252 0   119.95 0.276 0   119.95 0.3 0   119.95 0.324 0   119.95 0.348 0
              119.95 0.372 0   119.95 0.396 0   119.95 0.42 0   119.95 0.444 0   119.95 0.468 0
              119.95 0.492 0   119.95 0.516 0   119.95 0.54 0   119.95 0.564 0   119.95 0.588 0
              119.95 0.612 0   119.95 0.636 0   119.95 0.66 0   119.95 0.684 0   119.95 0.708 0
              119.95 0.732 0   119.95 0.756 0   119.95 0.78 0   119.95 0.804 0   119.95 0.828 0
              119.95 0.92 0
              '
    sort_by = y
  []
[]

[Executioner]
  type = SIMPLE

  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'

  momentum_l_abs_tol = 1e-14
  pressure_l_abs_tol = 1e-14
  turbulence_l_abs_tol = 1e-14
  momentum_l_tol = 1e-14
  pressure_l_tol = 1e-14
  turbulence_l_tol = 1e-14

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.4 0.4'
  turbulence_field_relaxation = '0.4 0.4'
  num_iterations = 5000
  pressure_absolute_tolerance = 1e-7
  momentum_absolute_tolerance = 1e-7
  turbulence_absolute_tolerance = '1e-7 1e-7'

  momentum_petsc_options_iname = '-u_system_pc_type -u_system_pc_hypre_type -v_system_pc_type -v_system_pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-TKE_system_pc_type -TKE_system_pc_hypre_type -TKED_system_pc_type -TKED_system_pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg hypre boomeramg'

  momentum_l_max_its = 300
  pressure_l_max_its = 300
  turbulence_l_max_its = 30

  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = true

  [csv]
    type = CSV
    execute_on = 'FINAL'
  []
[]
