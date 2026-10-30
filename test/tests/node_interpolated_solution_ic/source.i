# Writes a quadratic field on a TET10 mesh whose edges are curved by a nonlinear map of the nodes,
# including the mid-edge nodes.
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
    elem_type = TET10
  []
  [curve]
    type = ParsedNodeTransformGenerator
    input = gmg
    x_function = 'x + 0.08 * y * y'
    y_function = 'y + 0.1 * x * x'
    z_function = 'z + 0.05 * x * y'
  []
[]
[Problem]
  solve = false
[]
[AuxVariables]
  [T]
    order = SECOND
    [InitialCondition]
      type = FunctionIC
      function = '300 + 10 * x * y + 5 * z * z - 3 * x'
    []
  []
[]
[Executioner]
  type = Steady
[]
[Outputs]
  exodus = true
[]
