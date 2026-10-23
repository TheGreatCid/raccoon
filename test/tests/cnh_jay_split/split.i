# One HEX8 element under pure dilatation (each face moved outward by `rate * t`, negative for
# compression), with CNHJay elasticity and LargeDeformationJ2PlasticityBeBar. The deformation has no
# deviatoric part, so the isochoric energy W is zero and the elastic energy is the volumetric
# U(J) = (K/4)[(J - 1)^2 + (1/J - 1)^2]. With the energy split, U is active and degraded in tension
# (J > 1) and inactive and undegraded in compression (J < 1). The degradation is a constant g = 0.5,
# so the degraded and undegraded energies differ.

E = 200e3
nu = 0.3
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
g = 0.5
rate = -0.01
split = true

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    elem_type = HEX8
  []
[]

[Variables]
  [disp_x][]
  [disp_y][]
  [disp_z][]
[]

[AuxVariables]
  [d][]
[]

[Kernels]
  [x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
  [z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
  []
[]

[BCs]
  [x0]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [y0]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [z0]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [x1]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '${rate} * t'
  []
  [y1]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = '${rate} * t'
  []
  [z1]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = '${rate} * t'
  []
[]

[Materials]
  [defgrad]
    type = ComputeDeformationGradient
  []
  [constants]
    type = ADGenericConstantMaterial
    prop_names = 'K G g dg/dd'
    prop_values = '${K} ${G} ${g} 0'
  []
  [nodeg]
    type = NoDegradation
    phase_field = d
    property_name = nodeg
  []
  [elasticity]
    type = CNHJay
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
  []
  [J2]
    type = LargeDeformationJ2PlasticityBeBar
    phase_field = d
    hardening_model = hardening
    apply_strain_energy_split = ${split}
    # The hardening model's own energy names, so they do not clash with the elastic ones
    hardening_plastic_energy_density_active = psip_active
  []
  [hardening]
    type = LinearHardening
    yield_stress = 1e9
    hardening_modulus = 0
    phase_field = d
    degradation_function = nodeg
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = elasticity
    plasticity_model = J2
  []
[]

[Functions]
  # Volumetric energy of the imposed deformation, J = (1 + rate * t)^3
  [U]
    type = ParsedFunction
    expression = 'J := (1 + r * t)^3; K / 4 * ((J - 1)^2 + (1 / J - 1)^2)'
    symbol_names = 'r K'
    symbol_values = '${rate} ${K}'
  []
  [p_over_J]
    type = ParsedFunction
    expression = 'J := (1 + r * t)^3; K / 2 * (J - 1 / J) * (J + 1 / J - 1) / J'
    symbol_names = 'r K'
    symbol_values = '${rate} ${K}'
  []
[]
[Materials]
  [stress_xx]
    type = ADRankTwoCartesianComponent
    rank_two_tensor = stress
    property_name = stress_xx
    index_i = 0
    index_j = 0
  []
[]

[Postprocessors]
  [psie_active]
    type = ADElementAverageMaterialProperty
    mat_prop = psie_active
  []
  [psie]
    type = ADElementAverageMaterialProperty
    mat_prop = psie
  []
  [U]
    type = FunctionValuePostprocessor
    function = U
  []
  # Cauchy pressure and the volumetric Kirchhoff pressure over J it should equal, undegraded in
  # compression and degraded by g in tension when the split is on
  [pressure]
    type = ADElementAverageMaterialProperty
    mat_prop = stress_xx
  []
  [p_over_J]
    type = FunctionValuePostprocessor
    function = p_over_J
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  dt = 1
  end_time = 2
[]

[Outputs]
  csv = true
[]
