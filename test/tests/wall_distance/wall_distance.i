[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 3
    dx = '1 1 1'
    dy = '1 1 1'
    dz = '1 1 1'
    ix = '5 5 5'
    iy = '5 5 5'
    iz = '5 5 5'
    subdomain_id = '0 0 0
                    0 0 0
                    0 0 0

                    0 0 0
                    0 1 0
                    0 0 0

                    0 0 0
                    0 0 0
                    0 0 0'
  []
  [remove]
    type = BlockDeletionGenerator
    input = cmg
    block = 1
    new_boundary = 'inner'
  []
[]

[AuxVariables]
  [wall_distance]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [wall_distance_aux]
    type = WallDistanceAux
    variable = wall_distance
    walls = 'inner left right top bottom'
    execute_on = 'initial'
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
