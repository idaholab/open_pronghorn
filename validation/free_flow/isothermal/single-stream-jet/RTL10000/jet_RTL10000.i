# ================================================================================================================
# Model description
# ================================================================================================================
# Open-Pronghorn simulation of a steady state incompressible turbulent flow of a jet into a big space/tank.
# Reference: https://pubs.aip.org/aip/pof/article/36/11/115133/3319195/Calibration-of-the-Reynolds-stress-model-for
# Author: Dr. Kyriakopoulos Vasileios, Dr. Mauricio Tano
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
Vin =  '${fparse Re * mu / (rho * Dh)}' # Inlet velocity (m/s) # Uniform velocity profile Faghani et al. 2011
intensity = 0.05 # Turutoglu et al 2024
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
# The front boundary is the solid inlet-plane area outside the jet opening;
# the remaining outer tank boundaries are pressure outlets.
walls = 'front'
wall_treatment = 'neq' # Options: eq_newton, eq_incremental, eq_linearized, neq
eps = 1e-4
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
    preserve_volumes = on
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

[Physics]
  [NavierStokes]
    [FlowSegregated/flow]
      velocity_variable = 'vel_x vel_y vel_z'
      pressure_variable = 'pressure'

      initial_velocity = '0 0 0'
      initial_pressure = '1e-8'

      density = ${rho}
      dynamic_viscosity = ${mu}

      inlet_boundaries = 'inlet'
      momentum_inlet_types = 'fixed-velocity'
      momentum_inlet_functors = '0 0 ${Vin}'

      wall_boundaries = ${walls}
      momentum_wall_types = 'noslip'

      outlet_boundaries = ${outlet}
      momentum_outlet_types = 'fixed-pressure fixed-pressure fixed-pressure fixed-pressure fixed-pressure'
      pressure_functors = '0 0 0 0 0'

      include_symmetrized_viscous_stress = true
      orthogonality_correction = true
      pressure_two_term_bc_expansion = false
      momentum_two_term_bc_expansion = false
      momentum_advection_interpolation = ${advected_interp_method}
      system_names = 'u_system v_system w_system pressure_system'
      # consistency with the kernel syntax input, allows to use the same executioner block
      rhie_chow_uo_name = 'rc'
      pressure_projection_method = 'consistent'
    []
    [AdvancedKEpsilon/k-epsilon]
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
      C_mu = ${C_mu}

      # wall treatment
      turbulence_walls = ${walls}
      wall_treatment_eps = ${wall_treatment}
      bulk_wall_treatment = ${bulk_wall_treatment}

      # k-epsilon variant options
      k_epsilon_variant   = ${k_epsilon_variant}
      two_layer_flavor    = ${two_layer_flavor}
      use_buoyancy        = ${use_buoyancy}
      use_compressibility = ${use_compressibility}
      nonlinear_model     = ${nonlinear_model}
      curvature_model     = ${curvature_model}
      use_yap             = ${use_yap}
      use_low_re_Gprime   = ${use_low_re_Gprime}
      Ca0 = 0.667
      Ca1 = 1.25
      Ca2 = 1.0
      Ca3 = 0.9
    []
  []
[]

[FVInterpolationMethods]
  [upwind]
    type = FVAdvectedUpwind
  []
[]

[LinearFVBCs]
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
[]

[AuxVariables]
  [yplus]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
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
    u = vel_x
    v = vel_y
    w = vel_z
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []
[]

# The effective viscosities for the momentum, k and epsilon equations are created by the
# flow and turbulence Physics

[VectorPostprocessors]
  [centerline]
    type = LineValueSampler
    start_point = '${eps} ${eps} ${fparse 0.5 * L / nz}'
    end_point   = '${eps} ${eps} ${fparse L - 0.5 * L / nz}'
    num_points  = ${nz}
    variable    = vel_z
    sort_by     = 'z'
    execute_on  = 'FINAL'
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
  num_iterations = 1500
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
