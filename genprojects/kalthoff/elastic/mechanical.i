Gc = 22.2
l = 0.35
psic = 7.9
E = 1.9e5
nu = 0.3
rho = 8e-9
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
[GlobalParams]
  displacements = 'disp_x disp_y'
[]
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
    file = '../gold/kal_nd.msh'
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
  [psie_active]
    order = CONSTANT
    family = MONOMIAL
  []
  [d]
  []
[]
[AuxKernels]
  [stress]
    type = ADRankTwoScalarAux
    variable = stress
    rank_two_tensor = 'stress'
    scalar_type = MaxPrincipal
    execute_on = 'TIMESTEP_END'
  []
  [psie_active]
    type = ADMaterialRealAux
    variable = psie_active
    property = 'psie_active'
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
    density = reg_density
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
    function = 'if(t<1e-6, 0.5*4.00e10*t*t, 4.00e4*t-0.5*4.00e-2)'
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
  [elasticity]
    type = SmallDeformationIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d
    degradation_function = g
    decomposition = SPECTRAL
  []
  [strain]
    type = ADComputeSmallStrain
  []
  [stress]
    type = ComputeSmallDeformationStress
    elasticity_model = elasticity
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
  print_linear_residuals = false
  exodus = true
  interval = 1
[]
