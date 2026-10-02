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

# --- operating point ---

# --- corrosion operating point ---

# --- k-epsilon closure ---

# Cell-local EPI3V chemistry stage on the MSFR cavity mesh.
# Temperature is cell-local from outer T_fluid; dose is the original MSFR cosine spatial field.

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
  kernel_coverage_check = false
[]

[Variables]
  [e_sol]
    type = MooseVariableFVReal
  []
  [Cl_ion]
    type = MooseVariableFVReal
  []
  [Cl_rad]
    type = MooseVariableFVReal
  []
  [Cl2m_rad]
    type = MooseVariableFVReal
  []
  [Cl3_ion]
    type = MooseVariableFVReal
  []
  [Cl2_diss]
    type = MooseVariableFVReal
  []
  [Cr_II]
    type = MooseVariableFVReal
  []
  [Cr_III]
    type = MooseVariableFVReal
  []
  [Cr_I]
    type = MooseVariableFVReal
  []
  [T_local]
    type = MooseVariableFVReal
    initial_condition = 673.15
  []
  [dose_local]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
[]
[AuxVariables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 0.0
  []
  [pressure]
    type = INSFVPressureVariable
    initial_condition = 0.0
  []

  [D_eff_half]
    type = MooseVariableFVReal
    initial_condition = 1.0e-9
  []
  [a_u]
    type = MooseVariableFVReal
    initial_condition = 2.0
  []
  [a_v]
    type = MooseVariableFVReal
    initial_condition = 2.0
  []
[]
[MultiApps]
  [posttransport]
    type = TransientMultiApp
    sub_cycling = true
    input_files = epi3v_opt_msfr_cavity_dt0p00625_posttransport.i
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [species_to_posttransport]
    type = MultiAppCopyTransfer
    to_multi_app = posttransport
    source_variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
    variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
  []
  [velocity_to_posttransport]
    type = MultiAppCopyTransfer
    to_multi_app = posttransport
    source_variable = 'u v'
    variable = 'u v'
  []
  [species_from_posttransport]
    type = MultiAppCopyTransfer
    from_multi_app = posttransport
    source_variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
    variable = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
  []

  [D_eff_half_to_posttransport]
    type = MultiAppCopyTransfer
    to_multi_app = posttransport
    source_variable = D_eff_half
    variable = D_eff_half
  []

  [T_local_to_posttransport]
    type = MultiAppCopyTransfer
    to_multi_app = posttransport
    source_variable = T_local
    variable = T_local
  []
  [rc_fields_to_posttransport]
    type = MultiAppCopyTransfer
    to_multi_app = posttransport
    source_variable = 'pressure a_u a_v'
    variable = 'pressure a_u a_v'
  []
[]

[Executioner]
  type = Transient
  start_time = 0.0
  end_time = 0.2

  [TimeIntegrators]
    [epi3v]
      type = MSREPI3VTimeIntegrator
      salt_type = chloride
      metals = 'Cr'
      temperature = 673.15
      temperature_variable = T_local
      g_value_species = 'e_sol Cl_rad Cl2m_rad'
      g_value_overrides = '0.0 0.0 0.30'
      dose_rate = 5.0e5  # fallback only; local dose_local overrides this
      dose_rate_variable = dose_local
      species = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Cr_II Cr_III Cr_I'
      rtol = 1e-5
      atol = 1.0e-12
    []
  []

  [TimeStepper]
    type = ConstantDT
    dt = 0.00625
  []
[]

[Outputs]
  exodus = false
[]
