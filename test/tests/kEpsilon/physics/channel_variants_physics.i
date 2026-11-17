##########################################################
# ERCOFTAC test case for turbulent channel flow
# Case Number: 032
# Turbulent model using:
# k-epsilon, set up with the KEpsilonVariantsTurbulencePhysics
# Equilibrium + Newton wall treatment
# SIMPLE solve
##########################################################

### Problem Parameters ###
H = 1 # half-width of the channel
L = 120
Re = 2000
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
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = ${H}
    nx = 10
    ny = 5
    bias_y = 0.9
  []
  [block_2]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = ${fparse -H}
    ymax = 0
    nx = 10
    ny = 5
    bias_y = ${fparse 1/0.9}
  []
  [smg]
    type = StitchMeshGenerator
    inputs = 'block_1 block_2'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'bottom top'
    merge_boundaries_with_same_name = true
  []
  # Prevent test diffing on distributed parallel element numbering
  allow_renumbering = false
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

      initial_velocity = '${bulk_u} 0 0'
      initial_pressure = '1e-8'

      density = ${rho}
      dynamic_viscosity = ${mu}

      inlet_boundaries = 'left'
      momentum_inlet_types = 'fixed-velocity'
      momentum_inlet_functors = '${bulk_u} 0'

      wall_boundaries = ${walls}
      momentum_wall_types = 'noslip noslip'

      outlet_boundaries = 'right'
      momentum_outlet_types = 'fixed-pressure'
      pressure_functors = '0'

      include_symmetrized_viscous_stress = true
      orthogonality_correction = false
      pressure_two_term_bc_expansion = false
      momentum_two_term_bc_expansion = false
      momentum_advection_interpolation = ${advected_interp_method}
      system_names = 'u_system v_system pressure_system'
      rhie_chow_uo_name = 'rc'
      pressure_projection_method = 'consistent'
    []
    [AdvancedKEpsilon/k-epsilon]
      turbulence_handling = 'k-epsilon'
      tke_name = TKE
      tked_name = TKED
      system_names = 'TKE_system TKED_system'

      initial_tke = ${k_init}
      initial_tked = ${eps_init}

      mu_t_ratio_max = 1e20
      sigma_k = ${sigma_k}
      sigma_eps = ${sigma_eps}
      C_pl = 1e10
      C1_eps = ${C1_eps}
      C2_eps = ${C2_eps}
      C_mu = ${C_mu}

      turbulence_walls = ${walls}
      wall_treatment_eps = ${wall_treatment}
      bulk_wall_treatment = ${bulk_wall_treatment}
      use_nonorthogonal_correction = false

      # The k-epsilon variant is selected by the test specification. The wall distance needed by
      # the two-layer and low-Reynolds formulations is created by the Physics
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
  turbulence_equation_relaxation = '0.2 0.2'
  turbulence_field_relaxation = '0.2 0.2'
  num_iterations = 1000
  pressure_absolute_tolerance = 1e-7
  momentum_absolute_tolerance = 1e-7
  turbulence_absolute_tolerance = '1e-7 1e-7'

  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'

  momentum_l_max_its = 300
  pressure_l_max_its = 300
  turbulence_l_max_its = 30

  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = false
  [csv]
    type = CSV
    execute_on = FINAL
  []
[]

variables_to_sample = 'vel_x vel_y pressure TKE TKED'

[VectorPostprocessors]
  [line_center_channel]
    type = LineValueSampler
    start_point = '${fparse 0.125 * L} ${fparse 0.0001} 0'
    end_point = '${fparse 0.875 * L} ${fparse 0.0001} 0'
    num_points = ${Mesh/block_1/nx}
    variable = ${variables_to_sample}
    sort_by = 'x'
    execute_on = 'timestep_end'
  []
  [line_quarter_radius_channel]
    type = LineValueSampler
    start_point = '${fparse 0.125 * L} ${fparse 0.5 * H} 0'
    end_point = '${fparse 0.875 * L} ${fparse 0.5 * H} 0'
    num_points = ${Mesh/block_1/nx}
    variable = ${variables_to_sample}
    sort_by = 'x'
    execute_on = 'timestep_end'
  []
[]
