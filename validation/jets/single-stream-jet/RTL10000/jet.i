# ================================================================================================================
# Model description
# ================================================================================================================
# Open-Pronghorn simulation of a steady state incompressible turbulent flow of a jet into a big space/tank.
# Reference: https://pubs.aip.org/aip/pof/article/36/11/115133/3319195/Calibration-of-the-Reynolds-stress-model-for
# Author: Dr. Kyriakopoulos Vasileios
# Idaho Falls, INL, April 10, 2026
# ================================================================================================================
# MODEL PARAMETERS
# ================================================================================================================
# Geometry
# ================================================================================================================
H = 0.2 # Height/width of the tank (m).
L = 0.3 # Total tank length (m)
Dh = 0.01 # Inlet Pipe diameter (m)
nz = 150 # number of cells in the axial direction
n1 = 8 # number of radial cells per ring in the cross-section
n2 = 8 # number of radial cells per ring in the cross-section
n3 = 8 # number of radial cells per ring in the cross-section
n4 = 8 # number of radial cells per ring in the cross-section
n5 = 8 # number of radial cells per ring in the cross-section
ny = ${fparse n1 + n2 + n3 + n4 + n5}
# A = '${fparse pi*Dh*Dh/4}' # Cross-section area of the inlet (m2)
# ================================================================================================================
# Air properties
# ================================================================================================================
rho = 1.184 # Density (kg/m3)
# ================================================================================================================
# Flow and Boundary conditions
# ================================================================================================================
Re = 10000.0 # Reynolds number
mu = 1.844e-5 # Pa*s
Vin =  '${fparse Re * mu / (rho * Dh)}' # Inlet velocity (m/s) # Uniform velocity profile Faghani et al.
intensity = 0.05 # Turutoglu et al
k_init = '${fparse 1.5*(intensity * Vin)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (2*0.07*Dh)}'
advected_interp_method = 'upwind'
# ================================================================================================================
# k-epsilon Closure Parameters
# ================================================================================================================
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09
walls = 'front'
wall_treatment = 'neq' # Options: eq_newton, eq_incremental, eq_linearized, neq
# ================================================================================================================
# Modeling parameters: Turbulence-model knobs (optional)
# ================================================================================================================
k_epsilon_variant   = 'RealizableTwoLayer'  # Standard | StandardLowRe | StandardTwoLayer | Realizable | RealizableTwoLayer
two_layer_flavor    = 'Wolfstein' # Wolfstein | NorrisReynolds | Xu (only used for *TwoLayer variants)
use_buoyancy        = false
use_compressibility = false
nonlinear_model     = 'none'
curvature_model     = 'none'
use_yap             = false
use_low_re_Gprime   = false
bulk_wall_treatment = false
outlet = 'back bottom top left right'

[Mesh]
  [front_face]
    type = ConcentricCircleMeshGenerator
    num_sectors = 8
    radii = '${fparse Dh/2.0} ${Dh} ${fparse Dh*2.0} ${fparse Dh*4.0}'
    rings = '${n1} ${n2} ${n3} ${n4} ${n5}'
    has_outer_square = on
    pitch = ${H}
    preserve_volumes = off
    smoothing_max_it = 3
  []
  [extrude]
    type = AdvancedExtruderGenerator
    direction = '0 0 1'
    input = front_face
    heights = '${L}'
    num_layers = ${nz}
  []
  [inlet]
    type = ParsedGenerateSideset
    input = extrude
    combinatorial_geometry = 'abs(z) < 1e-6'
    included_subdomains = '1'
    normal = '0 0 -1'
    new_sideset_name = inlet
  []
  [front]
    type = ParsedGenerateSideset
    input = inlet
    combinatorial_geometry = 'abs(z) < 1e-6'
    included_subdomains = '2 3 4 5'
    normal = '0 0 -1'
    new_sideset_name = front
  []
  [back]
    type = ParsedGenerateSideset
    input = front
    combinatorial_geometry = 'abs(z) > ${fparse L - 1e-6}'
    included_subdomains = '1 2 3 4 5'
    normal = '0 0 1'
    new_sideset_name = back
  []
  [delete_surfaces]
    type = BoundaryDeletionGenerator
    input = back
    operation = remove
    boundary_names = '5 6'
  []
  [rename]
    type = RenameBlockGenerator
    input = delete_surfaces
    old_block = '1 2 3 4 5'
    new_block = 'air air air air air'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system TKE_system TKED_system'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    w = vel_z
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    pressure_projection_method = 'consistent'
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    initial_condition = 0
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    initial_condition = 0
    solver_sys = v_system
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    initial_condition = 0
    solver_sys = w_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    initial_condition = 1e-8
    solver_sys = pressure_system
  []
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
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'x'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = no
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = '${mu}'
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
  []

  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'y'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = no
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = '${mu}'
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
  []

  [w_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_z
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'z'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = no
  []
  [w_diffusion]
    type = LinearFVDiffusion
    variable = vel_z
    diffusion_coeff = '${mu}'
  []
  [w_pressure]
    type = LinearFVMomentumPressure
    variable = vel_z
    pressure = pressure
    momentum_component = 'z'
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

  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink]
    type = kEpsilonTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    w = vel_z
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10

    # NEW (optional)
    k_epsilon_variant        = ${k_epsilon_variant}    # 'Standard', 'Realizable', etc.
    use_buoyancy             = ${use_buoyancy}
    use_compressibility      = ${use_compressibility}
    nonlinear_model          = ${nonlinear_model}
    curvature_model          = ${curvature_model}
    Pr_t                     = 0.9
    C_M                      = 1.0
    gravity                  = '0 0 0'                 # leave 0 for BFS

    # if/when you have these fields:
    # temperature       = T
    # beta              = beta
    # speed_of_sound    = c
    # nonlinear_production = Gnl
    # curvature_factor  = fc
  []

  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_source_sink]
    type = kEpsilonTKEDSourceSink
    variable = TKED

    u = vel_x
    v = vel_y
    w = vel_z
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    # C_pl = 1e10

    # NEW (optional)
    k_epsilon_variant   = ${k_epsilon_variant}
    use_buoyancy        = ${use_buoyancy}
    use_compressibility = ${use_compressibility}
    nonlinear_model     = ${nonlinear_model}
    curvature_model     = ${curvature_model}
    use_yap             = ${use_yap}
    use_low_re_Gprime   = ${use_low_re_Gprime}

    Pr_t = 0.9
    C_M  = 1.0
    gravity = '0 -9.81 0'

    # same functors as for TKE if you use them:
    # temperature       = T
    # beta              = beta
    # speed_of_sound    = c
    # nonlinear_production = Gnl
    # curvature_factor  = fc
    wall_distance     = wall_distance   # for low-Re / two-layer Yap / G' terms
  []
[]

[LinearFVBCs]
  [inlet-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_x
    functor = '0.0'
  []
  [inlet-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_y
    functor = '0.0'
  []
  [inlet-w]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = vel_z
    functor = ${Vin}
  []
  [walls-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_x
    functor = 0.0
  []
  [walls-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_y
    functor = 0.0
  []
  [walls-w]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_z
    functor = 0.0
  []
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = ${outlet}
    variable = vel_x
    use_two_term_expansion = false
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = ${outlet}
    variable = vel_y
    use_two_term_expansion = false
  []
  [outlet_w]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = ${outlet}
    variable = vel_z
    use_two_term_expansion = false
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${outlet}
    variable = pressure
    functor = 0.0
  []

  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKE
    functor = '${k_init}'
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = ${outlet}
    variable = TKE
    use_two_term_expansion = false
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKED
    functor = '${eps_init}'
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = ${outlet}
    variable = TKED
    use_two_term_expansion = false
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = ${walls}
    variable = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
[]

[AuxVariables]
  [wall_distance]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
  []
  [yplus]
    type = MooseVariableFVReal
        two_term_boundary_expansion = false
  []
  [mu_eff]
    type = MooseVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
    two_term_boundary_expansion = false
  []
[]

[AuxKernels]
  [compute_wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = ${walls}
    execute_on = 'INITIAL NONLINEAR'
  []
  [compute_mu_t]
    type = kEpsilonViscosity
    variable = mu_t

    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    w = vel_z

    bulk_wall_treatment = ${bulk_wall_treatment}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e20
    execute_on = 'NONLINEAR'

    # NEW (optional) – choose model and options
    k_epsilon_variant = ${k_epsilon_variant}    # e.g. 'Standard' or 'Realizable'
    two_layer_flavor  = ${two_layer_flavor}     # ignored unless *TwoLayer variants
    Cd0 = 0.091      # defaults, can omit if you keep Standard
    Cd1 = 0.0042
    Cd2 = 0.00011
    Ca0 = 0.667      # Realizable C_mu coefficients
    Ca1 = 1.25
    Ca2 = 1.0
    Ca3 = 0.9
    wall_distance = wall_distance   # only needed for LowRe/TwoLayer (see below)
  []
  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    tke = TKE
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    w = vel_z
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []
  [compute_mu_eff]
    type = ParsedAux
    variable = 'mu_eff'
    coupled_variables = 'mu_t'
    expression = 'mu_t + ${mu}'
  []
[]

[VectorPostprocessors]
  [centerline]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '0 0 0.3'
    num_points = ${fparse nz + 1}
    variable = 'mu_eff mu_t pressure TKE TKED vel_x vel_y vel_z yplus'
    sort_by = 'z'
    execute_on = 'FINAL'
  []
  !include radial.i
[]

[Executioner]
  type = SIMPLE

  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'

  momentum_l_abs_tol = 1e-6
  pressure_l_abs_tol = 1e-6
  turbulence_l_abs_tol = 1e-6
  momentum_l_tol = 1e-6
  pressure_l_tol = 1e-6
  turbulence_l_tol = 1e-6

  momentum_equation_relaxation = 0.6
  pressure_variable_relaxation = 0.25
  turbulence_equation_relaxation = '0.4 0.4'
  turbulence_field_relaxation = '0.4 0.4'
  num_iterations = 10
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
  turbulence_l_max_its = 300

  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = true
  [out]
    type = CSV
    execute_on = FINAL
  []
[]
