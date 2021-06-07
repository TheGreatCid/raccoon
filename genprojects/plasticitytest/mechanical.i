a = 1
R = 4

psic = 9.5e5
#psic = 40e11
Gc = 1.38e8
l = '${fparse R * a}'


E = 68.8e9
nu = 0.33
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
beta = 0.5
sigma_y = 320e6
n = 5
ep0 = 0.01
v = '${fparse 0.5 * a}'

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
  [to_sub_coalescence_mobility]
    type = MultiAppCopyTransfer
    multi_app = fracture
    direction = to_multiapp
    source_variable = coalescence_mobility
    variable = coalescence_mobility
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    xmax = ${a}
    ymax = ${a}
    zmax = ${a}
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
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
  [solid_z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
  []
[]
[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []
 # [ydisp]
 #   type = FunctionDirichletBC
 #   variable = disp_y
 #   function = '1*t'
 #   boundary = top
 #  []

 [ydisp]
   type = LoadingUnloadingDirichletBC
   variable = disp_y
   boundary = 'top'
   initial_load_cap = '${fparse 0.005*a}'
   load_cap_increment = '0.005'
   load_step = '${fparse 0.0001 * a}'
   ultimate_load = 0.06
   unloaded_indicator = stress
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
    parameter_values = '2 1 0 1e-06'
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
  [hencky]
    type = HenckyIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    #decomposition = SPECTRAL
    output_properties = 'psie_active'
    outputs = exodus
  []
  [nodeg]
    type = NoDegradation
    phase_field = d
    f_name = nodeg
  []
  [power_law_hardening]
    type = PowerLawHardening
    yield_stress = ${sigma_y}
    exponent = ${n}
    reference_plastic_strain = ${ep0}
    phase_field = d
    degradation_function = nodeg
    output_properties = 'psip_active'
    outputs = exodus
  []
  [J2]
    type = LargeDeformationJ2Plasticity
    hardening_model = power_law_hardening
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
    plasticity_model = J2
  []

  [coalescence]
    type = ADParsedMaterial
    f_name = coalescence_mobility #mobility
    material_property_names = 'effective_plastic_strain'
    constant_names = 'beta ep0'
    constant_expressions = '${beta} ${ep0}'
    function = 1-(1-beta)*(1-exp(-(effective_plastic_strain/ep0)))
    #function = 1
    outputs = exodus
    output_properties = 'coalescence_mobility'
  []

[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true

line_search = none

  dt = '${fparse 0.0001 * a}'
  end_time = '${v}'

  abort_on_solve_fail = true
[]

[Outputs]
  file_base = stress_deformation
  print_linear_residuals = false
  csv = true
  exodus = true
[]
[Postprocessors]
  [F]
    type = ElementAverageValue
    variable = F
  []
  [stress]
    type = ElementAverageValue
    variable = stress
  []
  [d]
    type = ElementAverageValue
    variable = d
  []
  [ep]
    type = ADElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
  []
  [psie]
    type = ADElementAverageMaterialProperty
    mat_prop = psie_active
  []
  [psip]
    type = ADElementAverageMaterialProperty
    mat_prop = psip_active
  []
  [coal]
    type = ADElementAverageMaterialProperty
    mat_prop = coalescence_mobility
  []
[]
