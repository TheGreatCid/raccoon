# Recovers source.e (a cut crack in a stretched damage band) through SolutionUserObjectQP, which
# relaxes the damaged elements behind the crack faces: they start stress-free (stretch and dumped
# stress reset) and keep their mass. Reports the relaxed volume, the recovered stretch and stress
# ranges, and the total mass, which must not depend on the relaxation.

relax_threshold = 0.5

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = source.e
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  element = TET10_4th
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [disp_x]
    order = SECOND
  []
  [disp_y]
    order = SECOND
  []
  [disp_z]
    order = SECOND
  []
  [F_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [relaxed]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [F_yy]
    type = ADRankTwoAux
    variable = F_yy
    rank_two_tensor = deformation_gradient
    index_i = 1
    index_j = 1
    execute_on = INITIAL
  []
  [stress_yy]
    type = ADRankTwoAux
    variable = stress_yy
    rank_two_tensor = stress_sol
    index_i = 1
    index_j = 1
    execute_on = INITIAL
  []
  # The dumped stress is 1 everywhere, so a zero stress marks a relaxed element
  [relaxed]
    type = ParsedAux
    variable = relaxed
    coupled_variables = stress_yy
    expression = '1 - stress_yy'
    execute_on = INITIAL
  []
[]

[UserObjects]
  [epsol]
    type = SolutionUserObjectQP
    mesh = source.e
    system_variables = 'd'
    tensor_materials = 'rotation_tensor stretch_tensor stress'
    materials = ''
    timestep = LATEST
    nodal_variable_order = SECOND
    relax_threshold = ${relax_threshold}
    execute_on = INITIAL
  []
[]

[Materials]
  [defgrad]
    type = ComputeDeformationGradient
    recover = true
    solution = epsol
    recover_mode = polar_decomposition
  []
  [stress_sol]
    type = SolutionTensor
    solution = epsol
    tensor_name = stress
    relaxed_value = zero
  []
  [rho]
    type = ADGenericConstantMaterial
    prop_names = density
    prop_values = 1
  []
  [adj_density]
    type = ADStrainAdjustedDensityCustom
    strain_free_density = density
    base_name = adj
  []
[]

[Postprocessors]
  [relaxed_volume]
    type = ElementIntegralVariablePostprocessor
    variable = relaxed
    execute_on = INITIAL
  []
  [F_yy_max]
    type = ElementExtremeValue
    variable = F_yy
    execute_on = INITIAL
  []
  [F_yy_min]
    type = ElementExtremeValue
    variable = F_yy
    value_type = min
    execute_on = INITIAL
  []
  [stress_yy_min]
    type = ElementExtremeValue
    variable = stress_yy
    value_type = min
    execute_on = INITIAL
  []
  # Reference mass: sum of rho / J over the recovered state
  [mass]
    type = ADElementIntegralMaterialProperty
    mat_prop = adj_density
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
  csv = true
  [exodus]
    type = Exodus
    execute_on = INITIAL
  []
[]
