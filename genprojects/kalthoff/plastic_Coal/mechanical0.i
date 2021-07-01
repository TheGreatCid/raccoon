Gc = 22.2
l = 0.35
psic = 7.9
E = 1.9e5
nu = 0.3
rho = 8e-9
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
eta = 1
#sigma_y = 2000 #Check if this value makes sense
sigma_y = 350
n = 1 #for power law
ep0 = 0.2
beta = 0.2

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
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    alpha = '${fparse -1/3}'
    component = 0
  []
  [solid_y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    alpha = '${fparse -1/3}'
    component = 1
  []
[]

[BCs]
  [xdisp]
    type = FunctionDirichletBC
    variable = 'disp_x'
    boundary = 'load'
     #function = 'if(t<1e-7, 0.5*2.00e8*t*t, 20.0*t-1.00e-6)'
    #function = 'if(t<1e-7, 0.5*2.00e8*t*t, 20.0*t-1.00e-6)'
    #function = 'if(t<1e-8, 0.5*2.00e10*t*t, 2.0e4*t-1.00e-2)'

    #function = 'if(t<1e-6, 0.5*1.65e10*t*t, 1.65e4*t-0.5*1.65e-2)'
    function = 'if(t<1e-6, 0.5*3.20e10*t*t, 3.20e4*t-0.5*3.20e-2)'
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
  [defgrad]
    type = ComputeDeformationGradient
  []
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
    #function = 1
    outputs = exodus
    output_properties = 'coalescence_mobility'
  []
  [nodeg]
    type = NoDegradation
    phase_field = d
    f_name = nodeg
  []
  [degradation]
    type = RationalDegradationFunction
    f_name = g
    phase_field = d
    parameter_names = 'p a2 a3 eta'
    parameter_values = '2 1 0 1e-06'
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
   [hencky]
     type = HenckyIsotropicElasticity
     bulk_modulus = K
     shear_modulus = G
     phase_field = d
     degradation_function = g
     decomposition = SPECTRAL
     output_properties = 'elastic_strain psie_active'
     outputs = exodus
   []
  # [strain] #For elasticity
  #  type = ADComputeSmallStrain
  # []
  # [elasticity]
  #   type = SmallDeformationIsotropicElasticity
  #   bulk_modulus = K
  #   shear_modulus = G
  #   phase_field = d
  #   degradation_function = g
  #   decomposition = SPECTRAL
  #   output_properties = 'elastic_strain psie_active'
  #   outputs = exodus
  # []
  [J2]
    #type = SmallDeformationJ2Plasticity
    type = LargeDeformationJ2Plasticity
    hardening_model = power_law_hardening
    output_properties = 'effective_plastic_strain'
    outputs = exodus
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
  [stress]
    type = ComputeLargeDeformationStress
    #type = ComputeSmallDeformationStress
    elasticity_model = hencky
    plasticity_model = J2
  []
[]

[Executioner]
  type = Transient
  dt = 5e-7
  #end_time = 9e-5
  #dt = 5e-9
  end_time = 10.25e-5
  [TimeIntegrator]
    type = NewmarkBeta
    gamma = '${fparse 5/6}'
    beta = '${fparse 4/9}'
    # solve_type = lumped
    # use_constant_mass = true
  []

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  automatic_scaling = true
  # [Quadrature]
  #   order = CONSTANT
  # []
[]
[Outputs]
 file_base = 'exodusfiles/kalthoff/kal_plastic_v350_b02e02_oldprops_HHT'
  print_linear_residuals = false
  exodus = true
  interval = 1
[]
# [Postprocessors]
#
#   [disp_x]
#     type = AverageNodalVariableValue
#     boundary = 'load'
#     variable = disp_x
#
#   []
# []
