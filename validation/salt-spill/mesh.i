# Mesh Generation for Salt Spill Molten Pool Radiative experiment
dx = ${fparse 1/4*2.54/100}
int_rad = ${fparse 2*2.54/100 / 2}
ext_rad = ${fparse int_rad + dx}
height_salt = ${fparse 12.5/1000}
height_air = ${fparse 3*2.54/100 - height_salt}

[Mesh]
  [base]
    type = ConcentricCircleMeshGenerator
    num_sectors = 20
    radii = '${fparse 0.95*int_rad} ${int_rad} ${ext_rad}'
    rings = '5 5 8'
    has_outer_square = off
    pitch = 1.42063
    #portion = left_half
    preserve_volumes = off
    smoothing_max_it = 3
  []
  [extrusion]
    type = AdvancedExtruderGenerator
    input = base
    heights = '${dx} ${height_salt} ${height_air}'
    num_layers = '8 15 30'
    direction = '0 0 1'
  []
  [bot_lid]
    type = ParsedSubdomainMeshGenerator
    input = extrusion
    combinatorial_geometry = 'z < ${fparse dx}'
    block_id = '4'
    block_name = 'solid'
  []
  [solid]
    type = RenameBlockGenerator
    input = bot_lid
    old_block = '3'
    new_block = 'solid'
  []
  [salt]
    type = ParsedSubdomainMeshGenerator
    input = solid
    combinatorial_geometry = 'z < ${fparse dx+height_salt}'
    block_id = '3'
    block_name = 'salt'
    excluded_subdomains = 'solid'
  []
  [air]
    type = RenameBlockGenerator
    input = salt
    old_block = '1 2'
    new_block = 'air air'
  []
  # Boundaries
  [rename]
    type = RenameBoundaryGenerator
    input = air
    old_boundary = 'outer 2 3'
    new_boundary = 'external_wall solid_bottom top'
  []
  [solid_salt]
    type = SideSetsBetweenSubdomainsGenerator
    input = rename
    new_boundary = 'salt_solid_wall'
    primary_block = 'salt'
    paired_block = 'solid'
  []
  [salt_air]
    type = SideSetsBetweenSubdomainsGenerator
    input = solid_salt
    new_boundary = 'salt_air_iface'
    primary_block = 'salt'
    paired_block = 'air'
  []
  [solid_air]
    type = SideSetsBetweenSubdomainsGenerator
    input = salt_air
    new_boundary = 'air_solid_iface'
    primary_block = 'solid'
    paired_block = 'air'
  []
  [top_solid]
    type = ParsedGenerateSideset
    input = solid_air
    combinatorial_geometry = 'z > ${height_salt}'
    new_sideset_name = 'solid_top'
    included_boundaries = 'top'
    included_subdomains = 'solid'
  []
  [top_air]
    type = ParsedGenerateSideset
    input = top_solid
    combinatorial_geometry = 'z > ${height_salt}'
    new_sideset_name = 'air_top'
    included_boundaries = 'top'
    included_subdomains = 'air'
  []
[]

[Outputs]
  exodus = true
[]
