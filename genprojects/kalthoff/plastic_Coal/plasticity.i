Gc = 22.2
l = 0.35
psic = 7.9
E = 1.9e5
nu = 0.3
rho = 8e-9
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
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = true
[]

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = '../gold/kal.msh'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
[]

[AuxVariables]
  [stress]
    order = CONSTANT
    family = MONOMIAL
  []
  [d]
  []
  [effective_plastic_strain]
    order = CONSTANT
    family = MONOMIAL
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
  [inertia_x]
    type = InertialForce
    variable = disp_x
    density = 'reg_density'
  []
  [inertia_y]
    type = InertialForce
    variable = disp_y
    density = 'reg_density'
  []
  [solid_x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
  []
[]

[BCs]
  [xdisp]
    type = FunctionDirichletBC
    variable = 'disp_x'
    boundary = 'load'
    function = 'if(t<1e-6, 0.5*1.65e10*t*t, 1.65e4*t-0.5*1.65e-2)'
    preset = false
  []
  [y_bot]
    type = DirichletBC
    variable = 'disp_y'
    boundary = 'bottom'
    value = '0'
  []
[]

[Materials]
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc psic density'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic} ${rho}'
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
  [reg_density]
    type = MaterialConverter
    ad_props_in = 'density'
    reg_props_out = 'reg_density'
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
 file_base = 'kal_plasticity_coal'
  print_linear_residuals = false
  exodus = true
  interval = 1
[]
