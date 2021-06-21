
#Units from Miehe paper
Gc = 115e3
l = 2e-3
rho = 7800
lambda = 115.38e9
G = 76.93e9
nu = '${fparse lambda/(2*(lambda+G))}' #0.3
E = '${fparse 2*G*(1+nu)}' #2.00018e11
K = '${fparse E/3/(1-2*nu)}'

eta = 1
n = 1
#sigma_y = 2e9
sigma_y = 2e9
ep0 = 1
beta = 1
#psic = 3.03e6
psic = 6e6

#Thermal values
R = 8.3143 #Ideal gas constant
Q = 400e3 #Activation Energy, rough guess, 150e3 J/mole, another possible value
sigma0 = 2e9#Reference yield stress
kappa = 45e3 #W/k #45 W/mk
#W/mk = (J/S)/(mk)
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
    file = '../gold/kal_m.msh'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [temp]
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
  [temp_in_k]
  []
  [ref_temp]
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
  [temp_in_K]
      type = ParsedAux
      variable = 'temp_in_k'
      args = 'temp'
      function = 'temp + 273.15'
  []
  [reftemp]
    type = ParsedAux
    variable = 'ref_temp'
    function = '293.15'
  []
[]

[Kernels]
  [conduction]
    type = HeatConduction
    variable = 'temp'
  []
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
#ref temp


[BCs]
  [xdisp]
    type = FunctionDirichletBC
    variable = 'disp_x'
    boundary = 'load'
    function = 'if(t<1e-7, 0.5*2.00e8*t*t, 20.0*t-1.00e-6)'
    #function = 'if(t<1e-7, 0.5*1.65e8*t*t, 16.5*t-8.25e-7)'
    #function = 'if(t<1e-6, 0.5*1.65e10*t*t, 1.65e4*t-0.5*1.65e-2)'
    #function = 'if(t<1e-6, 0.5*2.00e10*t*t, 2.00e4*t-0.5*2.00e-2)'
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
    prop_names = 'K G l Gc psic density Q sigma0'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic} ${rho} ${Q} ${sigma0}'
  []
  [thermal_conduct]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '${kappa}'
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
  [eigenstrain]
    type = ComputeThermalExpansionEigenDeformationGradient
    reference_temperature = ref_temp
    temperature = temp
    eigen_deformation_gradient_name = thermal_defgrad
    thermal_expansion_function = '${fparse 11.2e-6}'
  []
  [defgrad]
    type = ComputeDeformationGradient
    eigen_deformation_gradient_names = 'thermal_defgrad'
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
  [hencky]
    type = HenckyIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    #decomposition = SPECTRAL
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
  []
  [J2]
    type = LargeDeformationJ2Plasticity
    hardening_model = arrhenius_law_hardening
    output_properties = 'effective_plastic_strain'
    outputs = exodus
  []
  [arrhenius_law]
    type = ArrheniusLaw
    arrhenius_coefficient = A
    activation_energy = ${Q}
    ideal_gas_constant = ${R}
    T = temp_in_k
  []
  [arrhenius_law_hardening]
    type = ArrheniusLawHardening
    reference_stress = sigma0
    arrhenius_coefficient = A
    eps = '${ep0}'
    phase_field = d
    degradation_function = nodeg
    output_properties = 'psip_active'
    outputs = exodus
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
    plasticity_model = J2
  []
[]

[Executioner]
  type = Transient
  #dt = 5e-7
  #end_time = 9e-5
  dt = 3e-8
  end_time = 10.25e-5
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
 file_base = 'exodusfiles/kalthoff/kal_thermoplastic_v200_e08_b08'
  print_linear_residuals = false
  exodus = true
  interval = 5
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
