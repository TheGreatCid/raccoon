# Builds a recovery file with a cut crack, as MOOSE2SIERRA's CUT/cut.i does before remeshing: a
# planar crack is cut into band.i's damage band and its faces pushed apart, the fields of band.e
# are copied, and the per-QP deformation (from band.i's crack-opening displacement, which
# stretches the band across the crack) and a constant stress are written as recovery variables.

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = band.e
  []
  [crack]
    type = PlanarCrackGenerator
    input = fmg
    solution_file = band.e
    threshold = 0.3
    lower_block_id = 10
    upper_block_id = 11
  []
  [cut]
    type = BreakMeshByBlockGenerator
    input = crack
    block_pairs = '10 11'
  []
  [free]
    type = RemoveDisjointNeighborPairsGenerator
    input = cut
  []
  [separate]
    type = SeparateCrackFacesGenerator
    input = free
    boundaries = 'Block10_Block11 Block11_Block10'
    separation = 1e-3
    moved_nodes_from = crack
  []
  [merge]
    type = RenameBlockGenerator
    input = separate
    old_block = '10 11'
    new_block = '0 0'
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  element = TET10_4th
[]

[Problem]
  solve = false
[]

[CopySolutionFields]
  file = band.e
  moved_nodes_from = separate
[]

[Materials]
  [defgrad]
    type = ComputeDeformationGradient
  []
  [stress]
    type = ADGenericConstantRankTwoTensor
    tensor_name = stress
    tensor_values = '0 0 0 0 1 0 0 0 0'
  []
[]

[RecoverVariables]
  [rec]
    tensor_materials = 'rotation_tensor stretch_tensor stress'
    materials = ''
    execute_on = INITIAL
  []
[]

[Executioner]
  type = Steady
  [Quadrature]
    order = FOURTH
  []
[]

[Outputs]
  [exodusqp]
    type = Exodus
    file_base = source
  []
[]
