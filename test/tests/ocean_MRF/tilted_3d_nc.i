### Thermophysical Properties ###
rho_0 = 1.
T_0 = 300.5
cp = 1.
beta = 0.1020408 #${fparse 1/9.8} # g.beta = 1
T_c = 300
T_h = 301

# Ra = 7.5e4
# Pr = 0.71
mu = 0.003076794869 #${fparse sqrt(Pr/Ra)}
k_cond = 0.0043335139 #${fparse sqrt(1/(Ra*Pr))}

[GlobalParams]
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  advected_interp_method = 'average'
  u = vel_x
  v = vel_y
  w = vel_z
[]

[Mesh]
  # final_generator = cube

  # [q1_low]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.0
  #   xmax = 0.5
  #   ymin = 0.0
  #   ymax = 0.5
  #   zmin = 0.0
  #   zmax = 0.5
  #   bias_x = 1.03
  #   bias_y = 1.03
  #   bias_z = 1.03
  # []

  # [q2_low]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.5
  #   xmax = 1.0
  #   ymin = 0.0
  #   ymax = 0.5
  #   zmin = 0.0
  #   zmax = 0.5
  #   bias_x = 0.97
  #   bias_y = 1.03
  #   bias_z = 1.03
  # []

  # [q3_low]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.0
  #   xmax = 0.5
  #   ymin = 0.5
  #   ymax = 1.0
  #   zmin = 0.0
  #   zmax = 0.5
  #   bias_x = 1.03
  #   bias_y = 0.97
  #   bias_z = 1.03
  # []

  # [q4_low]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.5
  #   xmax = 1.0
  #   ymin = 0.5
  #   ymax = 1.0
  #   zmin = 0.0
  #   zmax = 0.5
  #   bias_x = 0.97
  #   bias_y = 0.97
  #   bias_z = 1.03
  # []

  # [q1_high]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.0
  #   xmax = 0.5
  #   ymin = 0.0
  #   ymax = 0.5
  #   zmin = 0.5
  #   zmax = 1.0
  #   bias_x = 1.03
  #   bias_y = 1.03
  #   bias_z = 0.97
  # []

  # [q2_high]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.5
  #   xmax = 1.0
  #   ymin = 0.0
  #   ymax = 0.5
  #   zmin = 0.5
  #   zmax = 1.0
  #   bias_x = 0.97
  #   bias_y = 1.03
  #   bias_z = 0.97
  # []

  # [q3_high]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.0
  #   xmax = 0.5
  #   ymin = 0.5
  #   ymax = 1.0
  #   zmin = 0.5
  #   zmax = 1.0
  #   bias_x = 1.03
  #   bias_y = 0.97
  #   bias_z = 0.97
  # []

  # [q4_high]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 42
  #   ny = 42
  #   nz = 42
  #   xmin = 0.5
  #   xmax = 1.0
  #   ymin = 0.5
  #   ymax = 1.0
  #   zmin = 0.5
  #   zmax = 1.0
  #   bias_x = 0.97
  #   bias_y = 0.97
  #   bias_z = 0.97
  # []

  # [low_row1]
  #   type = StitchMeshGenerator
  #   inputs = 'q1_low q2_low'
  #   stitch_boundaries_pairs = 'right left'
  #   clear_stitched_boundary_ids = true
  # []

  # [low_row2]
  #   type = StitchMeshGenerator
  #   inputs = 'q3_low q4_low'
  #   stitch_boundaries_pairs = 'right left'
  #   clear_stitched_boundary_ids = true
  # []

  # [high_row1]
  #   type = StitchMeshGenerator
  #   inputs = 'q1_high q2_high'
  #   stitch_boundaries_pairs = 'right left'
  #   clear_stitched_boundary_ids = true
  # []

  # [high_row2]
  #   type = StitchMeshGenerator
  #   inputs = 'q3_high q4_high'
  #   stitch_boundaries_pairs = 'right left'
  #   clear_stitched_boundary_ids = true
  # []

  # [low_layer]
  #   type = StitchMeshGenerator
  #   inputs = 'low_row1 low_row2'
  #   stitch_boundaries_pairs = 'top bottom'
  #   clear_stitched_boundary_ids = true
  # []

  # [high_layer]
  #   type = StitchMeshGenerator
  #   inputs = 'high_row1 high_row2'
  #   stitch_boundaries_pairs = 'top bottom'
  #   clear_stitched_boundary_ids = true
  # []
  # [cube]
  #   type = StitchMeshGenerator
  #   inputs = 'low_layer high_layer'
  #   stitch_boundaries_pairs = 'front back'
  #   clear_stitched_boundary_ids = true
  # []
  [file_mesh]
    type = FileMeshGenerator
    file = tilted_3d_nc_out_fine_v2.e
    use_for_exodus_restart = true
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system energy_system'
  previous_nl_solution_required = true
[]

[UserObjects]
  [ins_rhie_chow_interpolator]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    w = vel_z
    pressure = pressure
    rho = ${rho_0}
    p_diffusion_kernel = p_diffusion
    #body_force_kernel_names = "u_omega; v_omega"
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    # initial_condition = 0.0
    initial_from_file_var = vel_x
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    # initial_condition = 0
    initial_from_file_var = vel_y
    solver_sys = v_system
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    # initial_condition = 0
    initial_from_file_var = vel_z
    solver_sys = w_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    # initial_condition = 1e-8
    solver_sys = pressure_system
    initial_from_file_var = pressure
  []
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    # initial_condition = 300.5
    initial_from_file_var = T_fluid
  []
[]

[AuxVariables]
[]

[AuxKernels]
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
    gravity = '0 -9.8 0'
    rho = ${rho_0}
    ref_temperature = ${T_0}
    alpha_name = ${beta}
    momentum_component = 'x'
    pitch_angle = pitch_angle
    roll_angle = roll_angle
    yaw_angle = yaw_angle
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
    gravity = '0 -9.8 0'
    rho = ${rho_0}
    ref_temperature = ${T_0}
    alpha_name = ${beta}
    momentum_component = 'y'
    pitch_angle = pitch_angle
    roll_angle = roll_angle
    yaw_angle = yaw_angle
  []

  [w_time]
    type = LinearFVTimeDerivative
    variable = vel_z
    factor = ${fparse rho_0}
  []
  [w_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_z
    mu = ${mu}
    momentum_component = 'z'
    use_nonorthogonal_correction = false
  []
  [w_pressure]
    type = LinearFVMomentumPressure
    variable = vel_z
    pressure = pressure
    momentum_component = 'z'
  []
  # [w_buoyancy]
  #   type = LinearFVSRFMomentumBoussinesq
  #   variable = vel_z
  #   T_fluid = T_fluid
  #   gravity = '0 -9.8 0'
  #   rho = ${rho_0}
  #   ref_temperature = ${T_0}
  #   alpha_name = ${beta}
  #   momentum_component = 'z'
  #   pitch_angle = pitch_angle
  #   roll_angle = roll_angle
  #   yaw_angle = yaw_angle
  # []

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
  [u_noslip]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_x
    boundary = 'back front bottom top left right'
    functor = 0
  []
  [v_noslip]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_y
    boundary = 'back front bottom top left right'
    functor = 0
  []
  [w_noslip]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_z
    boundary = 'back front bottom top left right'
    functor = 0
  []

  [pressure]
    type = LinearFVPressureFluxBC
    boundary = 'back front bottom top left right'
    variable = pressure
    u = vel_x
    v = vel_y
    w = vel_z
    rho = ${rho_0}
    HbyA_flux = HbyA
    Ainv = Ainv
  []

  [T_cold]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = T_fluid
    boundary = 'right left'
    functor = ${T_c}
  []
  [T_hot]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = T_fluid
    boundary = 'bottom top'
    functor = ${T_h}
  []
  [T_adiab]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = T_fluid
    boundary = 'front back'
    functor = 0.0
  []
[]

[FunctorMaterials]
  [SRF_Functor_Material]
    type = LinearFVSRFFunctorMaterial
    mc_origin = '0. 0. 0'
    SRF_input_mode = 'fixed'
    pitch_angle_fixed = 0.0
    roll_angle_fixed = 0.0
    yaw_angle_fixed = 45.0
    pitch_omega_fixed = 0.0
    roll_omega_fixed = 0.0
    yaw_omega_fixed = 0.0
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
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  energy_system = 'energy_system'
  momentum_equation_relaxation = 0.65
  pressure_variable_relaxation = 0.25
  energy_equation_relaxation = 0.85
  num_iterations = 5 #8000
  dt = 0.02
  num_steps = 6000
  pressure_absolute_tolerance = 1e-8
  momentum_absolute_tolerance = 1e-8
  energy_absolute_tolerance = 1e-8
  print_fields = false
  momentum_l_max_its = 300

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.5 0.5 0.5'

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
  [out_fine]
    type = Exodus
    time_step_interval = 50
  []
[]
