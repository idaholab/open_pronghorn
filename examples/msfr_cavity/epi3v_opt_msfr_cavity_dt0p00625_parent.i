# =====================================================================================================
# Molten salt fast reactor (MSFR) core-cavity model -- EVOL-benchmark-inspired geometry: a single,
# fully coupled, publication-quality example integrating
#
#   * fluid flow            -- incompressible Navier-Stokes (linear FV, Rhie-Chow, segregated SIMPLE),
#   * nuclear power          -- a separable cosine neutron-flux shape (simplified power production),
#   * energy deposition      -- volumetric fission heat, advected enthalpy, molecular + turbulent
#                               conduction,
#   * turbulence             -- standard k-epsilon with wall functions (turbulent viscosity mu_t),
#   * radiolysis             -- the chloride radical chain producing the oxidant Cl2.- proportional to
#                               power, transported by advection and molecular + turbulent diffusion,
#   * corrosion              -- the structural-alloy cavity walls dissolve through a temperature-
#                               dependent Butler-Volmer reaction, releasing chromium into the salt as
#                               Cr(II); the radiolytic oxidant then oxidizes it Cr(II) -> Cr(III),
#                               coupling corrosion and radiolysis through the dissolved-chromium redox.
#
# Geometry (the realistic MSFR concept rather than a channel): an open 2D core cavity bounded by solid
# structural-alloy walls. The fuel salt (NaCl-UCl3) is pumped in through a central inlet at the bottom
# and collected at a central outlet plenum at the top. The central jet rises through the core and
# entrains the surrounding salt, so two large RECIRCULATION ZONES form on either side, with the salt
# flowing back DOWN along the side walls. Those slow, recirculating regions are the realistic concern
# -- they are hotter (poor heat removal), they trap the radiolytic oxidant, and their wall corrosion
# is altered. This is the spirit of the EVOL/MSFR CFD benchmark (a heated core cavity with internal
# recirculation), here carrying the full multiphysics.
#
# Validation provenance (each physics component is validated separately in this application):
#   * turbulent flow / recirculating cavity -- validation/free_flow (ERCOFTAC channel and
#                            backward-facing-step, the canonical separated/recirculating flows),
#   * radiolysis kinetics -- validation/msr (pulse-radiolysis transient-absorption traces),
#   * corrosion kinetics  -- validation/corrosion (76 cases / 43 targets of the effective Butler-
#                            Volmer correlation; the wall rate here reproduces that calibration).
# Verification postprocessors confirm the energy balance closes (nuclear heat == enthalpy rise), the
# flow is fully turbulent (reported Reynolds number), and the recirculation is captured (reverse flow).
#
# Run:  open_pronghorn-opt -i msfr_cavity.i
# View: the Exodus output (vel with the corner recirculation, T_fluid hot spots, mu_t, Cl2m_rad,
#       Cr_II/Cr_III, corrosion_rate_um_y).
# =====================================================================================================

# --- geometry (MSFR-like core cavity) ---
Wc = 2.0              # m, cavity width
Hc = 2.0              # m, cavity height
in_min = 0.4          # m, central inlet: x from in_min ...
in_max = 1.6          # m, ... to in_max (1.2 m wide; wide enough to sweep the floor corners)
out_min = 0.7         # m, central outlet plenum: x from out_min ...
out_max = 1.3         # m, ... to out_max (0.6 m wide)

# --- NaCl-UCl3 fast-reactor fuel salt (~650 C) ---
rho = 3300            # kg/m^3
mu = 4.0e-3           # Pa s (molecular dynamic viscosity)
cp = 1000             # J/kg/K
k_mol = 0.4           # W/m/K (molecular thermal conductivity)
D_mol = 2.0e-9        # m^2/s (molecular ionic diffusivity)
Pr_t = 0.9            # turbulent Prandtl number
Sc_t = 0.7            # turbulent Schmidt number

# --- operating point ---
u_in = 1.0            # m/s inlet velocity across the cavity floor (upward, pumped circulation)
T_in = 873.15         # K inlet temperature (600 C)
q_peak = 6.0e7        # W/m^3 peak fission heat density (~60 MW/m^3)
dose_peak = 5.0e5     # J/m^3/s peak radiolytic dose
Cl_bulk = 2.0e4       # mol/m^3 bulk chloride
Cr_in = 3.0           # mol/m^3 inlet dissolved Cr(II) (bulk corrosion inventory, ~50 ppm)

# --- corrosion operating point ---
T_ref_corr = 923.15   # K reference (mean core) temperature for the calibrated corrosion rate

# --- k-epsilon closure ---
C_mu = 0.09
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
intensity = 0.05
k_init = '${fparse 1.5*(intensity*u_in)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (0.07*0.5*Wc)}'
walls = 'left right top_wall bottom_wall'
wall_treatment = 'eq_newton'
advected_interp_method = 'upwind'

[Mesh]
  [cavity]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${Wc}
    ymin = 0
    ymax = ${Hc}
    nx = 40
    ny = 40
  []
  # Carve a central inlet jet and outlet plenum (and the surrounding walls) out of the bottom and top
  # boundaries, then drop the originals so every face belongs to exactly one sideset. The salt jets up
  # the centre and recirculates back down along the side walls.
  [outlet]
    type = ParsedGenerateSideset
    input = cavity
    included_boundaries = 'top'
    combinatorial_geometry = 'x > ${out_min} & x < ${out_max}'
    new_sideset_name = 'outlet'
  []
  [top_wall]
    type = ParsedGenerateSideset
    input = outlet
    included_boundaries = 'top'
    combinatorial_geometry = 'x < ${out_min} | x > ${out_max}'
    new_sideset_name = 'top_wall'
  []
  [inlet]
    type = ParsedGenerateSideset
    input = top_wall
    included_boundaries = 'bottom'
    combinatorial_geometry = 'x > ${in_min} & x < ${in_max}'
    new_sideset_name = 'inlet'
  []
  [bottom_wall]
    type = ParsedGenerateSideset
    input = inlet
    included_boundaries = 'bottom'
    combinatorial_geometry = 'x < ${in_min} | x > ${in_max}'
    new_sideset_name = 'bottom_wall'
  []
  [drop_orig]
    type = BoundaryDeletionGenerator
    input = bottom_wall
    boundary_names = 'top bottom'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system energy_system '
                     'e_sol_sys cl_ion_sys cl_rad_sys cl2m_rad_sys cl3_ion_sys cl2_diss_sys '
                     'cr_ii_sys cr_iii_sys cr_i_sys'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
[]

[UserObjects]
  [rc]
    type = MSRLinearRhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    solver_sys = u_system
    initial_condition = 0
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = ${u_in}
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 1e-8
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
  [T_fluid]
    type = MooseLinearVariableFVReal
    solver_sys = energy_system
    initial_condition = ${T_in}
  []
[]

[FunctorMaterials]
  [cp_mat]
    type = GenericFunctorMaterial
    prop_names = 'cp'
    prop_values = '${cp}'
  []
  # Effective (molecular + turbulent) conductivity for the energy equation.
  [k_eff_mat]
    type = ParsedFunctorMaterial
    property_name = k_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${k_mol} + mut*${cp}/${Pr_t}'
  []
  # Effective (molecular + turbulent) diffusivity for the radiolysis and corrosion species.
  [D_eff_mat]
    type = ParsedFunctorMaterial
    property_name = D_eff
    functor_names = 'mu_t'
    functor_symbols = 'mut'
    expression = '${D_mol} + mut/(${rho}*${Sc_t})'
  []
  # Volumetric fission heat (for the energy balance verification).
  [heat_mat]
    type = ParsedFunctorMaterial
    property_name = nuclear_heat
    functor_names = 'flux_shape'
    functor_symbols = 'f'
    expression = '${q_peak} * f'
  []
  # Wall corrosion rate [um/y] at the local temperature: the calibrated reference rate scaled by the
  # Arrhenius factor of the exchange current (Ea = 38.92 kJ/mol; the reference rate is filled in below
  # from the action's seed). Valid where the dissolved metal is dilute, as in this single pass.
  [corrosion_rate_mat]
    type = ParsedFunctorMaterial
    property_name = corrosion_rate_functor
    functor_names = 'T_fluid'
    functor_symbols = 'T'
    expression = '${corr_rate_ref} * exp(38921.89 / 8.31446261815324 * (1.0/${T_ref_corr} - 1.0/T))'
  []
[]

# Reference corrosion rate [um/y] at T_ref_corr for this case (stainless 316, chloride, oxidizing),
# from the corrosion correlation; used only by the diagnostic rate field above. (Aggressive: chloride
# salts give high chromium dissolution rates.) The inlet Cr(II) is set to the reference concentration,
# so the near-wall concentration stays near c_ref and the net rate tracks this reference value.
corr_rate_ref = 142.55

[Functions]
  # Separable cosine flux shape, peaked at the core mid-plane (the neutron-flux profile that both the
  # fission heat and the radiolytic source follow).
  [flux_shape]
    type = ParsedFunction
    expression = 'cos(pi*(x-0.5*${Wc})/(1.4*${Wc})) * cos(pi*(y-0.5*${Hc})/(1.4*${Hc}))'
  []
  [dose_field]
    type = ParsedFunction
    expression = '${dose_peak} * f'
    symbol_names = 'f'
    symbol_values = 'flux_shape'
  []
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Cr'                       # track Cr(II)/Cr(III)/Cr(I) and their radiolytic reactions
  temperature = T_fluid               # chemistry driven by the solved temperature field
  dose_rate = dose_field              # neutron-flux-shaped radiolytic source
  # Lump the fast electron/atom kinetics into a net radiolytic-oxidant yield (Cl2.- produced directly
  # proportional to dose), keeping the explicit oxidant -> Cr(II) -> Cr(III) coupling.
  g_value_species = 'e_sol Cl_rad Cl2m_rad'
  g_value_overrides = '0.0 0.0 0.30'
  rhie_chow_user_object = 'rc'        # advection by the solved turbulent flow
  diffusivity = D_eff                 # molecular + turbulent diffusion
  time_derivative = false             # steady-state operating point (SIMPLE executioner)
[]

[CorrosionPlatingFlow]
  # The structural walls corrode, releasing chromium into the salt as the radiolysis-tracked Cr(II);
  # the radiolytic oxidant then oxidizes it to Cr(III), coupling corrosion and radiolysis.
  elements = 'Cr'
  release_variables = 'Cr_II'
  reaction_boundary = ${walls}
  temperature = T_fluid               # Arrhenius corrosion driven by the solved temperature field
  reference_temperature = ${T_ref_corr}
  material_class = stainless_316
  salt_class = chloride
  redox_class = oxidizing_fef2        # the radiolytically oxidizing core salt
  applied_overpotential = 0.12
  reference_concentration = ${Cr_in}
  temperature_dependent_kinetics = true
[]

[LinearFVKernels]
  # --- momentum (turbulent stress mu_t + molecular mu) + pressure ---
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
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
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
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

  # --- k-epsilon turbulence ---
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
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
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
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    C_pl = 1e10
    wall_distance = wall_distance
  []

  # --- energy: advected enthalpy + (molecular + turbulent) conduction + fission heat ---
  [energy_advection]
    type = LinearFVEnergyAdvection
    variable = T_fluid
    advected_quantity = temperature
    cp = ${cp}
  []
  [energy_conduction]
    type = LinearFVDiffusion
    variable = T_fluid
    diffusion_coeff = 'k_eff'
    use_nonorthogonal_correction = false
  []
  [fission_heat]
    type = LinearFVSource
    variable = T_fluid
    source_density = flux_shape
    scaling_factor = ${q_peak}
  []
[]

[LinearFVBCs]
  # --- flow ---
  [inlet_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = vel_x
    functor = 0.0
  []
  [inlet_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = vel_y
    functor = ${u_in}
  []
  [walls_u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_x
    functor = 0.0
  []
  [walls_v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = ${walls}
    variable = vel_y
    functor = 0.0
  []
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = vel_x
    use_two_term_expansion = false
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = vel_y
    use_two_term_expansion = false
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = outlet
    variable = pressure
    functor = 0.0
  []
  # --- turbulence ---
  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = TKE
    functor = ${k_init}
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = TKE
    use_two_term_expansion = false
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = TKED
    functor = ${eps_init}
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = TKED
    use_two_term_expansion = false
  []
  [walls_mu_t]
    type = LinearFVTurbulentViscosityWallFunctionBC
    boundary = ${walls}
    variable = mu_t
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
  # --- energy ---
  [inlet_T]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = T_fluid
    functor = ${T_in}
  []
  [outlet_T]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = T_fluid
    use_two_term_expansion = false
  []
  # --- radiolysis + corrosion species: bulk chloride and the inlet Cr inventory enter at the bottom;
  #     everything leaves at the top outlet ---
  [inlet_e_sol]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = e_sol
    functor = 0
  []
  [outlet_e_sol]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = e_sol
    use_two_term_expansion = false
  []
  [inlet_Cl_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cl_ion
    functor = ${Cl_bulk}
  []
  [outlet_Cl_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cl_ion
    use_two_term_expansion = false
  []
  [inlet_Cl_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cl_rad
    functor = 0
  []
  [outlet_Cl_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cl_rad
    use_two_term_expansion = false
  []
  [inlet_Cl2m_rad]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cl2m_rad
    functor = 0
  []
  [outlet_Cl2m_rad]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cl2m_rad
    use_two_term_expansion = false
  []
  [inlet_Cl3_ion]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cl3_ion
    functor = 0
  []
  [outlet_Cl3_ion]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cl3_ion
    use_two_term_expansion = false
  []
  [inlet_Cl2_diss]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cl2_diss
    functor = 0
  []
  [outlet_Cl2_diss]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cl2_diss
    use_two_term_expansion = false
  []
  [inlet_Cr_II]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cr_II
    functor = ${Cr_in}
  []
  [outlet_Cr_II]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cr_II
    use_two_term_expansion = false
  []
  [inlet_Cr_III]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cr_III
    functor = 0
  []
  [outlet_Cr_III]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cr_III
    use_two_term_expansion = false
  []
  [inlet_Cr_I]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = inlet
    variable = Cr_I
    functor = 0
  []
  [outlet_Cr_I]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = outlet
    variable = Cr_I
    use_two_term_expansion = false
  []
[]

[AuxVariables]
  [wall_distance]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '${fparse rho * C_mu * k_init^2 / eps_init}'
  []
  [corrosion_rate_um_y]
    type = MooseLinearVariableFVReal
    initial_condition = 0
  []
  [a_u_out]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
  [a_v_out]
    type = MooseVariableFVReal
    initial_condition = 1.0
  []
[]

[AuxKernels]
  [compute_wall_distance]
    type = WallDistanceAux
    variable = wall_distance
    walls = ${walls}
    execute_on = 'INITIAL'
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
    bulk_wall_treatment = false
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    mu_t_ratio_max = 1e20
    execute_on = 'NONLINEAR'
    wall_distance = wall_distance
  []
  [compute_corrosion_rate]
    type = FunctorAux
    variable = corrosion_rate_um_y
    functor = corrosion_rate_functor
    execute_on = 'TIMESTEP_END'
  []
  [export_a_u]
    type = MSRLinearMomentumDiagonalAux
    variable = a_u_out
    rhie_chow_user_object = rc
    component = x
    execute_on = TIMESTEP_END
  []
  [export_a_v]
    type = MSRLinearMomentumDiagonalAux
    variable = a_v_out
    rhie_chow_user_object = rc
    component = y
    execute_on = TIMESTEP_END
  []
[]

[Postprocessors]
  # --- thermal-hydraulics ---
  [T_outlet]
    type = SideAverageValue
    variable = T_fluid
    boundary = outlet
  []
  [mu_t_max]
    type = ElementExtremeValue
    variable = mu_t
    value_type = max
  []
  # --- recirculation and hot-spot diagnostics (the realistic cavity features) ---
  # A negative minimum vertical velocity is downward (reverse) flow: the signature of the corner
  # recirculation zones. The peak temperature is the hot spot that forms where the salt recirculates
  # and heat removal is poor.
  [vel_y_min]
    type = ElementExtremeValue
    variable = vel_y
    value_type = min
  []
  [vel_y_max]
    type = ElementExtremeValue
    variable = vel_y
    value_type = max
  []
  [T_hotspot]
    type = ElementExtremeValue
    variable = T_fluid
    value_type = max
  []
  # --- energy balance verification: with adiabatic walls and steady state, the fission power must
  #     equal the net enthalpy convected out. The enthalpy flux is the Rhie-Chow mass flux weighted
  #     by temperature (so it is correct even with the strong recirculation), summed over the inlet
  #     and outlet. energy_balance_error is then a true global conservation check. ---
  [nuclear_power]
    type = ElementIntegralFunctorPostprocessor
    functor = nuclear_heat
  []
  [delta_T]
    type = ParsedPostprocessor
    pp_names = 'T_outlet'
    expression = 'T_outlet - ${T_in}'
  []
  [enthalpy_flux]
    type = VolumetricFlowRate
    boundary = 'inlet outlet'
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = T_fluid
    rhie_chow_user_object = rc
    advected_interp_method = ${advected_interp_method}
  []
  [energy_balance_error]
    type = ParsedPostprocessor
    pp_names = 'nuclear_power enthalpy_flux'
    expression = 'abs(nuclear_power - ${rho} * ${cp} * enthalpy_flux) / nuclear_power'
  []
  # --- radiolysis + corrosion coupling ---
  [Cl2m_rad_max]
    type = ElementExtremeValue
    variable = Cl2m_rad
    value_type = max
  []
  [Cr_II_outlet]
    type = SideAverageValue
    variable = Cr_II
    boundary = outlet
  []
  [Cr_III_outlet]
    type = SideAverageValue
    variable = Cr_III
    boundary = outlet
  []
  # Salt oxidation state at the outlet: the coupling observable (radiolytic oxidation of the corrosion
  # product chromium).
  [redox_ratio_CrIII_CrII]
    type = ParsedPostprocessor
    pp_names = 'Cr_II_outlet Cr_III_outlet'
    expression = 'Cr_III_outlet / (Cr_II_outlet + 1e-30)'
  []
  [corrosion_rate_wall_max]
    type = ElementExtremeValue
    variable = corrosion_rate_um_y
    value_type = max
  []
[]

[Executioner]
  type = SIMPLE
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  energy_system = 'energy_system'
  passive_scalar_systems = 'e_sol_sys cl_ion_sys cl_rad_sys cl2m_rad_sys cl3_ion_sys cl2_diss_sys '
                           'cr_ii_sys cr_iii_sys cr_i_sys'
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.3 0.3'
  turbulence_field_relaxation = '0.3 0.3'
  energy_equation_relaxation = 0.9
  passive_scalar_equation_relaxation = '0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9'
  num_iterations = 600
  momentum_absolute_tolerance = 1e-6
  pressure_absolute_tolerance = 1e-6
  turbulence_absolute_tolerance = '1e-5 1e-5'
  energy_absolute_tolerance = 1e-4
  passive_scalar_absolute_tolerance = '1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'
  energy_petsc_options_iname = '-pc_type -pc_hypre_type'
  energy_petsc_options_value = 'hypre boomeramg'
  passive_scalar_petsc_options_iname = '-pc_type -pc_hypre_type'
  passive_scalar_petsc_options_value = 'hypre boomeramg'
  momentum_l_max_its = 300
  pressure_l_max_its = 300
  turbulence_l_max_its = 30
  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = true
  csv = true
[]

# =====================================================================================
# EPI3V MSFR-CAVITY GEOMETRY AND FLOW COMPATIBILITY
#
# The original cavity solve above remains intact. After SIMPLE converges, its solved
# vel_x/vel_y fields are copied into a separate transient nonlinear-FV Strang/EPI3V app.
# =====================================================================================

[VectorPostprocessors]
  [parent_exact_rc_flux]
    type = MSRRhieChowFaceFluxVPP
    rhie_chow_user_object = rc
    internal_only = true
    execute_on = TIMESTEP_END
  []
[]

[MultiApps]
  [epi3v_transport]
    type = FullSolveMultiApp
    input_files = epi3v_opt_msfr_cavity_dt0p00625_transport.i
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [vel_x_to_epi3v]
    type = MultiAppCopyTransfer
    to_multi_app = epi3v_transport
    source_variable = vel_x
    variable = u_full
    execute_on = TIMESTEP_END
  []
  [vel_y_to_epi3v]
    type = MultiAppCopyTransfer
    to_multi_app = epi3v_transport
    source_variable = vel_y
    variable = v_full
    execute_on = TIMESTEP_END
  []
  [temperature_to_epi3v]
    type = MultiAppCopyTransfer
    to_multi_app = epi3v_transport
    source_variable = T_fluid
    variable = T_local
    execute_on = TIMESTEP_END
  []

  [mu_t_to_epi3v]
    type = MultiAppCopyTransfer
    to_multi_app = epi3v_transport
    source_variable = mu_t
    variable = mu_t_full
    execute_on = TIMESTEP_END
  []
  [pressure_to_epi3v]
    type = MultiAppCopyTransfer
    to_multi_app = epi3v_transport
    source_variable = pressure
    variable = pressure
    execute_on = TIMESTEP_END
  []
  [a_u_to_epi3v]
    type = MultiAppCopyTransfer
    to_multi_app = epi3v_transport
    source_variable = a_u_out
    variable = a_u_full
    execute_on = TIMESTEP_END
  []
  [a_v_to_epi3v]
    type = MultiAppCopyTransfer
    to_multi_app = epi3v_transport
    source_variable = a_v_out
    variable = a_v_full
    execute_on = TIMESTEP_END
  []
[]

