[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = 'gold/fields_refined.e'
    use_for_exodus_restart = true
  []
[]

[Variables]
  [d]
  []
[]

[AuxVariables]
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [psii_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [bounds_dummy]
  []
  [Gc]
    initial_from_file_var = 'Gc'
  []
  [psic]
    initial_from_file_var = 'psic'
  []
[]

[Bounds]
  [irreversibility]
    type = VariableOldValueBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
  [upper]
    type = ConstantBounds
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1
  []
[]

[Kernels]
  [pff_diff]
    type = ADPFFDiffusion
    variable = d
  []
  [pff_react]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]

[Materials]
  [fracture_toughness]
    type = ADParsedMaterial
    property_name = Gc
    coupled_variables = Gc
    expression = 'Gc'
  []
  [critial_fracture_energy]
    type = ADParsedMaterial
    property_name = psic
    coupled_variables = psic
    expression = 'psic'
  []
  [length_scale]
    type = ADGenericConstantMaterial
    prop_names = 'l'
    prop_values = '${l}'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    property_name = alpha
    expression = 'd'
    phase_field = d
  []
  [degradation]
    type = RationalDegradationFunction
    property_name = g
    expression = (1-d)^p/((1-d)^p+(Gc/psic*xi/c0/l)*d*(1+a2*d+a2*a3*d^2))*(1-eta)+eta
    phase_field = d
    parameter_names = 'p a2 a3 eta '
    parameter_values = '2 1 0 1e-6'
    material_property_names = 'Gc psic xi c0 l '
  []
  [psi]
    type = ADDerivativeParsedMaterial
    property_name = psi
    expression = 'alpha*Gc/c0/l+g*(psie_active+psii_active)'
    coupled_variables = 'd psie_active psii_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l'
    derivative_order = 1
  []
[]

[BCs]
  [Periodic]
    [all]
      variable = d
      auto_direction = 'x y'
    []
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package -snes_type'
  petsc_options_value = 'lu       superlu_dist                  vinewtonrsls'

  nl_abs_tol = 1e-08
  nl_rel_tol = 1e-06
[]

[Outputs]
  print_linear_residuals = false
[]
