Gc = 22.2
l = 0.35
psic = 7.9
E = 1.9e5
nu = 0.3
k = 2e-4
rho = 8e-9
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
lambda = '${fparse K-2*G/3}'


[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[MultiApps]
  [fracture]
    type = TransientMultiApp
    input_files = 'fracture.i'
    cli_args = 'Gc=${Gc};l=${l};k=${k};psic=${psic}'
    execute_on = 'TIMESTEP_END'
  []
[]

[Transfers]
  [from_d]
    type = MultiAppCopyTransfer
    multi_app = fracture
    direction = from_multiapp
    source_variable = d
    variable = d
  []
  [to_psie_active]
    type = MultiAppCopyTransfer
    multi_app = fracture
    direction = to_multiapp
    variable = psie_active
    source_variable = psie_active
  []
[]
[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = '../gold/half_notched_plate_63.msh'
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
    variable = 'disp_x'
    density = 'reg_density'
  []
  [inertia_y]
    type = InertialForce
    variable = 'disp_y'
    density = 'reg_density'
  []
  [solid_x]
    type = ADStressDivergenceTensors
    variable = 'disp_x'
    component = 0
    displacements = 'disp_x disp_y'
  []
  [solid_y]
    type = ADStressDivergenceTensors
    variable = 'disp_y'
    component = 1
    displacements = 'disp_x disp_y'
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
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = SPECTRAL
    output_properties = 'elastic_strain psie_active'
    outputs = exodus
  []
  [strain]
    type = ADComputeSmallStrain
    displacements = 'disp_x disp_y'
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
    output_properties = 'stress'
    outputs = exodus
  []
  [degradation]
    type = RationalDegradationFunction
    f_name = g
    phase_field = d
    parameter_names = 'p a2 a3 eta'
    parameter_values = '2 1 0 1e-09'
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
[]

[Executioner]
  type = Transient
  dt = 5e-9
  end_time = 9e-5

  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  file_base = 'kaltoff_elastic'
  exodus = true
  interval = 20
[]
