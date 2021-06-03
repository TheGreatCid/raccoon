[Mesh]
  [meshgen]
    type = GeneratedMeshGenerator
    dim = 1
    xmax = 1
    xmin = 0
    nx = 100
  []
[]

[Variables]
  [d]
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [psip_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [coalescence_mobility]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Bounds]
  [irreversibility]
    type = VariableOldValueBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = lower
  []
  [upper]
    type = ConstantBoundsAux
    variable = bounds_dummy
    bounded_variable = d
    bound_type = upper
    bound_value = 1
  []

[]

[Kernels]
  [diff]
    type = ADPFFDiffusion
    variable = d
    fracture_toughness = Gc
    regularization_length = l
    normalization_constant = c0

  []
  [source]
    type = ADPFFSource
    variable = d
    free_energy = psi
  []
[]


[Materials]
  [fracture_properties]
    type = ADGenericConstantMaterial
    prop_names = 'Gc psic l'
    prop_values = '${Gc} ${psic} ${l}'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    f_name = alpha
    function = 'd'
    phase_field = d
  []
  [degradation]
    type = RationalDegradationFunction
    f_name = g
    phase_field = d
    parameter_names = 'p a2 a3 eta'
    parameter_values = '2 1 0 1e-04'
  []
  [psi]
    type = ADDerivativeParsedMaterial
    f_name = psi
    function = 'alpha*Gc*coalescence_mobility/c0/l+g*(psie_active+psip_active)'
    args = 'd psie_active psip_active'
    material_property_names = 'alpha(d) g(d) Gc c0 l coalescence_mobility'
    derivative_order = 1
  []
  [mobility]
    type = ADParsedMaterial
    args = 'coalescence_mobility'
    function = 'coalescence_mobility'
    f_name = 'coalescence_mobility'
    constant_on = ELEMENT
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -snes_type'
  petsc_options_value = 'lu vinewtonrsls'
  nl_abs_tol = 1e-08
  nl_rel_tol = 1e-06
  automatic_scaling = true
  [Quadrature]
    order = CONSTANT
  []
[]
[Outputs]
  print_linear_residuals = false
[]
