E = 4000
nu = 0.25
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
lambda = '${fparse K-2*G/3}'
Gc = 1e-3
l = 0.03
k = 2e-4
psic = 0.0017578125
sigma_y = 3.8e10
n = 5
v = '${fparse -sqrt(Gc*3/lambda)}'
ep0 = 0.8
[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'fracture.i'
    cli_args = 'Gc=${Gc};l=${l};psic=${psic}'
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_d]
    type = MultiAppMeshFunctionTransfer
    multi_app = fracture
    direction = from_multiapp
    variable = d
    source_variable = d
  []
  [to_psie_active]
    type = MultiAppMeshFunctionTransfer
    multi_app = fracture
    direction = to_multiapp
    variable = psie_active
    source_variable = psie_active
  []
  [to_psip_active]
    type = MultiAppMeshFunctionTransfer
    multi_app = fracture
    direction = to_multiapp
    variable = psip_active
    source_variable = psip_active
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = '../gold/domainTriTens.msh'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [bounds_dummy]
  []
  [stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [F]
    order = CONSTANT
    family = MONOMIAL
  []
  [d]
  []
[]

[AuxKernels]
  [stress]
    type = ADRankTwoAux
    variable = stress
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
  []
  [F]
    type = ADRankTwoAux
    variable = F
    rank_two_tensor = deformation_gradient
    index_i = 1
    index_j = 1
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
  []
[]
[BCs]
  [forcing]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'Top'
    function = '${v}*t'
    preset = false
  []
  [FixedHole_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'Hole'
    value = 0
  []
  [FixedHole_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'Hole'
    value = 0
  []
[]
[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc psic'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic}'
  []
  [degradation]
    type = RationalDegradationFunction
    f_name = g
    phase_field = d
    parameter_names = 'p a2 a3 eta'
    parameter_values = '2 1 0 1e-04'
  []
  [crack_geometric]
    type = CrackGeometricFunction
    f_name = alpha
    function = 'd'
    phase_field = d
  []
  [defgrad]
    type = ComputeDeformationGradient
  []
  # [hencky]
  #   type = HenckyIsotropicElasticity
  #   bulk_modulus = K
  #   shear_modulus = G
  #   phase_field = d
  #   degradation_function = g
  #   #decomposition = SPECTRAL
  #   output_properties = 'psie_active'
  #   outputs = exodus
  # []
  [CNHIso]
    type = CNHIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = VOLDEV
    output_properties = 'psie_active'
    outputs = exodus
  []
  [power_law_hardening]
    type = PowerLawHardening
    yield_stress = ${sigma_y}
    exponent = ${n}
    reference_plastic_strain = ${ep0}
    phase_field = d
    degradation_function = g
    output_properties = 'psip_active'
    outputs = exodus
  []
  [J2]
    type = LargeDeformationJ2Plasticity
    hardening_model = power_law_hardening
    output_properties = 'effective_plastic_strain'
    outputs = exodus
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = CNHIso
    plasticity_model = J2
  []

[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'


  dt = 0.01
  end_time = 8
  nl_abs_tol = 1e-06
  nl_rel_tol = 1e-06
  automatic_scaling = true
  compute_scaling_once = false
  fp_max_its = 0
  fp_tol = 1e-06
  accept_on_max_fp_iteration = true
[]

[Outputs]
  file_base = ./exodusfiles/FiberMatrix/Fiber_comp_EPPD
  print_linear_residuals = false
  interval = 10
  csv = true
  exodus = true

[]
# [Postprocessors]
#   [F]
#     type = ElementAverageValue
#     variable = F
#   []
#   [stress]
#     type = ElementAverageValue
#     variable = stress
#   []
#   [d]
#     type = ElementAverageValue
#     variable = d
#   []
#   [ep]
#     type = ADElementAverageMaterialProperty
#     mat_prop = effective_plastic_strain
#   []
#   [psie]
#     type = ADElementAverageMaterialProperty
#     mat_prop = psie_active
#   []
#   [psip]
#     type = ADElementAverageMaterialProperty
#     mat_prop = psip_active
#   []
#   [coal]
#     type = ADElementAverageMaterialProperty
#     mat_prop = coalescence_mobility
#   []
# []
