# First in-MOOSE verification of the MSREPI3VTimeIntegrator.
#
# Scope:
#   - one-cell homogeneous (0D) chloride chemistry
#   - Zn electron capture
#   - no dose source, transport, diffusion, or gas exchange
#   - fixed MOOSE timestep dt = 5e-12 s
#
# Expected fully coupled endpoint at t = 2e-10 s:
#   e_sol ~ 0.437418
#   Zn_I  ~ 0.562582

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
    xmax = 1.0
  []
[]

[Problem]
  kernel_coverage_check = false
[]

[MoltenSaltRadiolysis]
  salt_type = chloride
  metals = 'Zn'
  integration_method = epi3v
  temperature = 673.15
  dose_rate = 0.0

  initial_condition_species = 'e_sol Zn_II Cl_ion'
  initial_condition_values = '1.0 100.0 5000.0'

  verbose = true
[]

[Postprocessors]
  [c_e_sol]
    type = ElementAverageValue
    variable = e_sol
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [c_Zn_I]
    type = ElementAverageValue
    variable = Zn_I
    execute_on = 'INITIAL TIMESTEP_END'
  []

  [c_Zn_II]
    type = ElementAverageValue
    variable = Zn_II
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  start_time = 0.0
  end_time = 2.0e-10
  dt = 5.0e-12

  [TimeIntegrators]
    [epi3v]
      type = MSREPI3VTimeIntegrator
      salt_type = chloride
      metals = 'Zn'
      temperature = 673.15
      species = 'e_sol Cl_ion Cl_rad Cl2m_rad Cl3_ion Cl2_diss Zn_II Zn_I'
    []
  []
[]

[Outputs]
  csv = true
[]
