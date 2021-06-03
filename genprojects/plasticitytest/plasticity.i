Gc = 1.38e5
l = 0.1
psic = 330e3
E = 68.8e6
nu = 0.3
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
eta = 1
sigma_y = 2000 #Check if this value makes sense
n = 1 #for power law
ep0 = 0.345
beta = .9


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
  displacements = 'disp_x'
  #volumetric_locking_correction = true
[]

[Mesh]
  [meshgen]
    type = GeneratedMeshGenerator
    dim = 1
    xmax = 1
    xmin = 0
    nx = 1
  []
[]

[Variables]
  [disp_x]
  []
[]

[AuxVariables]
  [stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xx]
  []
  [d]
  []
  [effective_plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
  [fy]
  []
[]

[AuxKernels]
  [stress]
    type = ADRankTwoScalarAux
    variable = 'stress'
    rank_two_tensor = 'stress'
    scalar_type = 'MaxPrincipal'
    execute_on = 'TIMESTEP_END'
  []
[]

[Kernels]
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    save_in = fy
  []
[]

[BCs]
  [forcing]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = '0.1*t'
    preset = false
  []
  [fixed_y]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
[]

[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc psic'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic}'
  []
  [coalescence]
    type = ADParsedMaterial
    f_name = coalescence_mobility #mobility
    material_property_names = 'effective_plastic_strain'
    constant_names = 'beta ep0'
    constant_expressions = '${beta} ${ep0}'
    function = 1-(1-beta)*(1-exp(-(effective_plastic_strain/ep0)))
    outputs = exodus
    output_properties = 'coalescence_mobility'
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

  [strain]
    type = ADComputeSmallStrain
  []
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = NONE
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
  []
  [plasticity]
    type = SmallDeformationJ2Plasticity
    hardening_model = power_law_hardening
    output_properties = 'effective_plastic_strain'
    outputs = exodus
  []
  [power_law_hardening]
    type = PowerLawHardening
    degradation_function = g
    yield_stress = ${sigma_y}
    exponent = ${n}
    reference_plastic_strain = ${ep0}
    phase_field = d
    output_properties = 'psip_active'
    outputs = exodus
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    plasticity_model = plasticity
    output_properties = 'stress'
    outputs = exodus
  []
[]

[Executioner]
  type = Transient
  dt = 5e-7
  end_time = 9e-5
  # [TimeIntegrator]
  #   type = CentralDifference
  #   solve_type = lumped
  #   use_constant_mass = true
  # []
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  [TimeIntegrator]
    type = NewmarkBeta
  []
  [Quadrature]
    order = CONSTANT
  []
[]
[Outputs]
 file_base = 'plastest'
  print_linear_residuals = false
  exodus = true
  interval = 1
[]
[Postprocessors]
  [fy]
    type = NodalSum
    variable = fy
    boundary = left
  []
  [stress]
    type = ElementAverageValue
    variable = stress_xx
  []
[]
