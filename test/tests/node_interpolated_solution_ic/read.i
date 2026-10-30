# Reads the field of source.i back on the same (curved) mesh with SolutionIC and with
# NodeInterpolatedSolutionIC. Where SolutionIC can locate every point, both give the same nodal
# values, including at the mid-edge nodes.
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = source_out.e
  []
[]
[Problem]
  solve = false
[]
[UserObjects]
  [sol]
    type = SolutionUserObject
    mesh = source_out.e
    system_variables = 'T'
    nodal_variable_order = SECOND
    timestep = LATEST
  []
[]
[AuxVariables]
  [T_solution_ic]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      solution_uo = sol
      from_variable = T
    []
  []
  [T_node_interpolated]
    order = SECOND
    [InitialCondition]
      type = NodeInterpolatedSolutionIC
      solution_uo = sol
      from_variable = T
    []
  []
  [difference]
    order = SECOND
  []
[]
[AuxKernels]
  [difference]
    type = ParsedAux
    variable = difference
    coupled_variables = 'T_solution_ic T_node_interpolated'
    expression = 'T_node_interpolated - T_solution_ic'
    execute_on = INITIAL
  []
[]
[Postprocessors]
  [max_difference]
    type = NodalExtremeValue
    variable = difference
    value_type = max_abs
  []
  [T_max]
    type = NodalExtremeValue
    variable = T_node_interpolated
  []
[]
[Executioner]
  type = Steady
[]
[Outputs]
  csv = true
[]
