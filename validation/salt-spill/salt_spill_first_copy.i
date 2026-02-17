#rho_salt_s=2000.
rho_salt_l=2000.
#cp_salt_s=2090.
cp_salt_l=2090.
#k_salt_s=0.6
k_salt_l=0.001
L = 425000
#alpha_b = 1.2e-4
T_solidus = 730.
T_liquidus = '${fparse T_solidus + 0.1}'

rho_solid=3000
cp_solid=600
k_solid=10.


[Problem]
  kernel_coverage_check = false
  linear_sys_names = 'energy_system u_system v_system w_system p_system'
[]

[GlobalParams]
  view_factor_object_name = rt_vf
[]

[Mesh]
 [file]
    type = FileMeshGenerator
    file = mesh_ramiro.e
  []
  [patch_bot]
    type = PatchSidesetGenerator
    boundary = 'air_solid_iface'
    n_patches = 2
    input = file
  []
  [patch_wall]
    type = PatchSidesetGenerator
    boundary = 'salt_air_iface'
    n_patches = 3
    input = patch_bot
    partitioner = grid
    #centroid_partitioner_direction = x
  []
  [top_deletion]
    type = BoundaryDeletionGenerator
    input = patch_wall
    boundary_names = 'air_solid_iface salt_air_iface'
  []
[]


salt_air_iface = 'salt_air_iface_0 salt_air_iface_1 salt_air_iface_2 salt_air_iface_3'# salt_air_iface_4'

# salt_air_iface_5 salt_air_iface_6
#                   salt_air_iface_7 salt_air_iface_8 salt_air_iface_9 salt_air_iface_10 salt_air_iface_11 salt_air_iface_12 salt_air_iface_13
#                   salt_air_iface_14 salt_air_iface_15 salt_air_iface_16 salt_air_iface_17 salt_air_iface_18 salt_air_iface_19 salt_air_iface_20
#                   salt_air_iface_21 salt_air_iface_22 salt_air_iface_23 salt_air_iface_24 salt_air_iface_25 salt_air_iface_26 salt_air_iface_27
#                   salt_air_iface_28 salt_air_iface_29
#                   salt_air_iface_30 salt_air_iface_31 salt_air_iface_32 salt_air_iface_33 salt_air_iface_34 salt_air_iface_35 salt_air_iface_36
#                   salt_air_iface_37 salt_air_iface_38 salt_air_iface_39 salt_air_iface_40 salt_air_iface_41 salt_air_iface_42 salt_air_iface_43
#                   salt_air_iface_44 salt_air_iface_45 salt_air_iface_46 salt_air_iface_47 salt_air_iface_48 salt_air_iface_49 salt_air_iface_50
#                   salt_air_iface_51 salt_air_iface_52 salt_air_iface_53 salt_air_iface_54 salt_air_iface_55 salt_air_iface_56 salt_air_iface_57
#                   salt_air_iface_58 salt_air_iface_59'
                  # salt_air_iface_60 salt_air_iface_61 salt_air_iface_62 salt_air_iface_63 salt_air_iface_64
                  # salt_air_iface_65 salt_air_iface_66 salt_air_iface_67 salt_air_iface_68 salt_air_iface_69 salt_air_iface_70 salt_air_iface_71
                  # salt_air_iface_72 salt_air_iface_73 salt_air_iface_74 salt_air_iface_75 salt_air_iface_76 salt_air_iface_77 salt_air_iface_78
                  # salt_air_iface_79 salt_air_iface_80 salt_air_iface_81 salt_air_iface_82 salt_air_iface_83 salt_air_iface_84 salt_air_iface_85
                  # salt_air_iface_86 salt_air_iface_87 salt_air_iface_88 salt_air_iface_89'

air_solid_iface = 'air_solid_iface_0 air_solid_iface_1'# air_solid_iface_2 air_solid_iface_3 air_solid_iface_4 air_solid_iface_5 air_solid_iface_6
                  # air_solid_iface_7 air_solid_iface_8 air_solid_iface_9'# air_solid_iface_10 air_solid_iface_11 air_solid_iface_12 air_solid_iface_13
                  #  air_solid_iface_14 air_solid_iface_15 air_solid_iface_16 air_solid_iface_17 air_solid_iface_18 air_solid_iface_19 air_solid_iface_20
                  #  air_solid_iface_21 air_solid_iface_22 air_solid_iface_23 air_solid_iface_24 air_solid_iface_25 air_solid_iface_26 air_solid_iface_27
                  #  air_solid_iface_28 air_solid_iface_29'

emissivities = '0.9 0.9 0.9 0.9 0.9 0.9 0.9'# 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9'# 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9'
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9
                # 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9'

[Variables]
  [temperature]
    type = MooseLinearVariableFVReal
    solver_sys = 'energy_system'
    block = 'solid salt'
  []
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_system'
    initial_condition = 0
    block = salt
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = 'v_system'
    initial_condition = 0
    block = salt
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    solver_sys = 'w_system'
    initial_condition = 0
    block = salt
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = 'p_system'
    initial_condition = 0
    block = salt
  []
[]

[AuxVariables]
  [fl]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = salt
  []
  # [density]
  #   type = MooseVariableFVReal
  # []
  # [th_cond]
  #   type = MooseVariableFVReal
  # []
  # [cp_var]
  #   type = MooseVariableFVReal
  # []
  # [darcy_coef]
  #   type = MooseVariableFVReal
  # []
  # [fch_coef]
  #   type = MooseVariableFVReal
  # []
[]

[AuxKernels]
  [compute_fl]
    type = NSLiquidFractionAux
    variable = fl
    temperature = temperature
    T_liquidus = '${T_liquidus}'
    T_solidus = '${T_solidus}'
    execute_on = 'TIMESTEP_END'
    block = salt
  []
  # [rho_out]
  #   type = FunctorAux
  #   functor = 'rho_mixture'
  #   variable = 'density'
  # []
  # [th_cond_out]
  #   type = FunctorAux
  #   functor = 'k_mixture'
  #   variable = 'th_cond'
  # []
  # [cp_out]
  #   type = FunctorAux
  #   functor = 'cp_mixture'
  #   variable = 'cp_var'
  # []
  # [darcy_out]
  #   type = FunctorAux
  #   functor = 'Darcy_coefficient'
  #   variable = 'darcy_coef'
  # []
  # [fch_out]
  #   type = FunctorAux
  #   functor = 'Forchheimer_coefficient'
  #   variable = 'fch_coef'
  # []
[]

[FVICs]
  [ic_u_1]
    type = FVConstantIC
    variable = temperature
    value = 301.15
    block = 'solid'
  []
  [ic_u_2]
    type = FVConstantIC
    variable = temperature
    value = 770
    block = 'salt'
  []
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    w = vel_z
    pressure = pressure
    rho = 1
    p_diffusion_kernel = p_diffusion
    block = salt
  []
[]

[LinearFVKernels]
  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    block = salt
    use_nonorthogonal_correction = false
  []

  [temp_time_salt]
    type = LinearFVTimeDerivative
    variable = temperature
    factor = ${fparse rho_salt_l*cp_salt_l}
    block = salt
  []
  [temp_conduction_salt]
    type = LinearFVDiffusion
    diffusion_coeff = ${k_salt_l}
    variable = temperature
    block = salt
  []
  [temp_phasechange_source]
    type = LinearFVPhaseChangeSource
    variable = temperature
    L = ${L}
    T_liquidus = ${T_liquidus}
    T_solidus = ${T_solidus}
    rho = ${rho_salt_l}
    block = salt
  []

  [temp_time_solid]
    type = LinearFVTimeDerivative
    variable = temperature
    factor = ${fparse rho_solid*cp_solid}
    block = solid
  []
  [temp_conduction_solid]
    type = LinearFVDiffusion
    diffusion_coeff = ${k_solid}
    variable = temperature
    block = solid
  []
[]

[UserObjects]
  [./gray_lambert]
    type = ViewFactorObjectSurfaceRadiation
    boundary = '${air_solid_iface} air_top ${salt_air_iface}'
    fixed_temperature_boundary = 'air_top'
    fixed_boundary_temperatures = '300'
    emissivity = ${emissivities}
    temperature = temperature
    execute_on = 'LINEAR TIMESTEP_BEGIN TIMESTEP_END NONLINEAR'
  [../]

  [vf_study]
    type = ViewFactorRayStudy
    execute_on = INITIAL
    boundary = '${air_solid_iface} air_top ${salt_air_iface}'
    face_order = CONSTANT
    polar_quad_order = 8
    azimuthal_quad_order = 10
    face_type = GAUSS
    warn_non_planar = false
  []

  [rt_vf]
    type = RayTracingViewFactor
    boundary = '${air_solid_iface} air_top ${salt_air_iface}'
    execute_on = INITIAL
    ray_study_name = vf_study
    normalize_view_factor = true
  []
[]

[RayBCs]
  [vf]
    type = ViewFactorRayBC
    boundary = '${air_solid_iface} air_top ${salt_air_iface}'
  []
[]

[LinearFVBCs]
  # temperature BCs
  [top_solid_temp]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = temperature
    boundary = 'solid_top'
    functor = 301.15
  []
  [outer_solid_temp]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = temperature
    boundary = 'external_wall'
    functor = 301.15
  []
  [./radiation_solidwalls]
    type = LinearFVGrayLambert
    variable = temperature
    temperature_radiation = temperature
    coeff_diffusion = ${k_solid}
    surface_radiation_object_name = gray_lambert
    boundary = ${air_solid_iface}
  [../]
  [./radiation_salt_air_iface]
    type = LinearFVGrayLambert
    variable = temperature
    temperature_radiation = temperature
    coeff_diffusion = ${k_salt_l}
    surface_radiation_object_name = gray_lambert
    boundary = ${salt_air_iface}
  [../]
[]

[FunctorMaterials]
  # [ins_fv]
  #   type = INSFVEnthalpyFunctorMaterial
  #   rho = rho_mixture
  #   cp = cp_mixture
  #   temperature = 'T'
  # []
  # [eff_cp]
  #   type = NSFVMixtureFunctorMaterial
  #   phase_2_names = '${cp_salt_s} ${k_salt_s} ${rho_salt_s}'
  #   phase_1_names = '${cp_salt_l} ${k_salt_l} ${rho_salt_l}}'
  #   prop_names = 'cp_mixture k_mixture rho_mixture'
  #   phase_1_fraction = fl
  # []
  # [mushy_zone_resistance]
  #   type = INSFVMushyPorousFrictionFunctorMaterial
  #   liquid_fraction = 'fl'
  #   mu = '${mu}'
  #   rho_l = '${rho_liquid}'
  #   dendrite_spacing_scaling = 1e-1
  # []
  # [darcy_coeff_friction]
  #   type = ParsedFunctorMaterial
  #   property_name = 'darcy_coef_friction'
  #   functor_names = 'darcy_coef'
  #   functor_symbols = 'darcy_coef'
  #   expression = 'darcy_coef'
  # []
  # [forch_coeff_friction]
  #   type = ParsedFunctorMaterial
  #   property_name = 'forch_coef_friction'
  #   functor_names = 'U fch_coef'
  #   functor_symbols = 'U fch_coef'
  #   expression = 'fch_coef * U'
  # []
[]

[Executioner]
  type = PIMPLE
  num_iterations = 20
  dt = 1
  end_time = 30
  should_solve_momentum = false
  should_solve_pressure = false
  energy_system = 'energy_system'
  energy_l_abs_tol = 5e-8
  energy_l_tol = 1e-10
  energy_equation_relaxation = 0.8
  energy_field_relaxation = 0.8
  energy_absolute_tolerance = 5e-8
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  print_fields = false
  continue_on_max_its = true

  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'p_system'
[]

[Outputs]
  exodus = true
[]

[Postprocessors]
  [air_top]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = gray_lambert
    return_type = RADIOSITY
    boundary = 'air_top'
  []
  [air_solid_iface_0]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = gray_lambert
    return_type = RADIOSITY
    boundary = 'air_solid_iface_0'
  []
  [air_solid_iface_1]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = gray_lambert
    return_type = RADIOSITY
    boundary = 'air_solid_iface_1'
  []
  [salt_air_iface_0]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = gray_lambert
    return_type = RADIOSITY
    boundary = 'salt_air_iface_0'
  []
  [salt_air_iface_1]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = gray_lambert
    return_type = RADIOSITY
    boundary = 'salt_air_iface_1'
  []
  [salt_air_iface_2]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = gray_lambert
    return_type = RADIOSITY
    boundary = 'salt_air_iface_2'
  []
  [salt_air_iface_3]
    type = GrayLambertSurfaceRadiationPP
    surface_radiation_object_name = gray_lambert
    return_type = RADIOSITY
    boundary = 'salt_air_iface_3'
  []
[]
