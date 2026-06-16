# Reference simulation for the HEX8 F-bar comparison recover/remesh test.
#
# Tuned to invoke heavy volumetric locking on a HEX8 mesh so that F-bar
# correction is doing real work and the recovered stretch / U-bar fields
# carry information that the restart leg has to reproduce.
#
# Tunings vs recover_remesh/reference.i (the baseline HEX8 recover test):
#   - Poisson ratio nu = 0.49           (near-incompressible -> high K/G ratio
#                                        -> volumetric locking severe)
#   - sigma_0 = 0.1 (was 1)             (lower yield -> earlier and stronger
#                                        plastic flow -> isochoric plastic
#                                        deformation drives det F variation
#                                        between QPs that F-bar must average)
#   - final_velocity = 0.5 (was 0.05)   (10x deformation magnitude)
#   - end_time = 2.5 (was 1.2)          (more steps to accumulate F-bar history)
#   - The output exports include `Fnobar` and `deformation_gradient`
#     (= F-bar) so the per-QP F vs F-bar gap is visible in Paraview.
#
# Otherwise mirrors recover_remesh/reference.i: HEX8, THIRD-order quadrature
# (8 QPs/elem), J2 plasticity with Johnson-Cook hardening, temperature, and
# a phase-field fracture sub-app.

E = 201.8e3
nu = 0.49
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 7900
Gc = 6
l = 0.2
psic = 15
Q = 0.9
specific_heat = 4.47e-4
thermal_conductivity = 4.4e-4
c1 = 0.1
c2 = 0.9
c3 = 0.1

Tinit = 293

hht_alpha = -0.25
newmark_beta = '${fparse (1-hht_alpha)^2/4}'
newmark_gamma = '${fparse 1/2-hht_alpha}'

end_time = 1.5
dt = 0.1

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
  element = HEX8_3rd
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 3
    ny = 3
    nz = 1
    xmin = 0
    xmax = 3
    ymin = 0
    ymax = 3
    zmin = 0
    zmax = 1
    elem_type = HEX8
  []
[]

# [MultiApps]
#   [fracture]
#     type = TransientMultiApp
#     input_files = 'fracture_ref.i'
#     cli_args = 'Gc=${Gc};l=${l};out_file=reference_hex8_fbar_d'
#     execute_on = 'TIMESTEP_END'
#     clone_parent_mesh = no
#   []
# []

# [Transfers]
#   [to_coal]
#     type = MultiAppCopyTransfer
#     variable = coal
#     source_variable = coal
#     to_multi_app = fracture
#   []
#   [from_d]
#     type = MultiAppCopyTransfer
#     variable = d
#     source_variable = d
#     from_multi_app = fracture
#   []
#   [to_psie]
#     type = MultiAppCopyTransfer
#     variable = psie_corr_active
#     source_variable = psie_corr_active
#     to_multi_app = fracture
#   []
#   [to_psip]
#     type = MultiAppCopyTransfer
#     variable = psip_active
#     source_variable = psip_active
#     to_multi_app = fracture
#   []
# []

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [T]
    initial_condition = ${Tinit}
  []
[]

[AuxVariables]
  [d]
  []
  [d_old]
  []
  [d_corr]
  []

  [accel_x]
  []
  [vel_x]
  []
  [accel_y]
  []
  [vel_y]
  []
  [accel_z]
  []
  [vel_z]
  []
[]

[AuxKernels]
  [accel_x]
    type = NewmarkAccelAux
    variable = accel_x
    displacement = disp_x
    velocity = vel_x
    beta = ${newmark_beta}
    execute_on = 'timestep_end'
  []
  [vel_x]
    type = NewmarkVelAux
    variable = vel_x
    acceleration = accel_x
    gamma = ${newmark_gamma}
    execute_on = 'timestep_end'
  []
  [accel_y]
    type = NewmarkAccelAux
    variable = accel_y
    displacement = disp_y
    velocity = vel_y
    beta = ${newmark_beta}
    execute_on = 'timestep_end'
  []
  [vel_y]
    type = NewmarkVelAux
    variable = vel_y
    acceleration = accel_y
    gamma = ${newmark_gamma}
    execute_on = 'timestep_end'
  []
  [accel_z]
    type = NewmarkAccelAux
    variable = accel_z
    displacement = disp_z
    velocity = vel_z
    beta = ${newmark_beta}
    execute_on = 'timestep_end'
  []
  [vel_z]
    type = NewmarkVelAux
    variable = vel_z
    acceleration = accel_z
    gamma = ${newmark_gamma}
    execute_on = 'timestep_end'
  []
[]

[Kernels]
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    density = density
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_x
    acceleration = accel_x
    absolute_value_vector_tags = 'ref'
  []
  [inertia_y]
    type = ADInertialForce
    variable = disp_y
    density = density
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_y
    acceleration = accel_y
    absolute_value_vector_tags = 'ref'
  []
  [inertia_z]
    type = ADInertialForce
    variable = disp_z
    density = density
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_z
    acceleration = accel_z
    absolute_value_vector_tags = 'ref'
  []
  [x]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    absolute_value_vector_tags = 'ref'
  []
  [y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    absolute_value_vector_tags = 'ref'
  []
  [z]
    type = ADDynamicStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    absolute_value_vector_tags = 'ref'
  []
  [hcond_time]
    type = ADHeatConductionTimeDerivative
    variable = T
    density_name = density
    specific_heat = specific_heat
    absolute_value_vector_tags = 'ref'
  []
  [hcond]
    type = ADHeatConduction
    variable = T
    thermal_conductivity = thermal_conductivity
    absolute_value_vector_tags = 'ref'
  []
  [heat_source]
    type = ADCoefMatSource
    variable = T
    coefficient = -1
    prop_names = 'plastic_heat_generation'
    absolute_value_vector_tags = 'ref'
  []
[]

[RecoverVariables]
  [rec]
    # stretch_tensor (= U raw) is written so the restart side can switch
    # between recovering U_bar directly vs recovering raw U + re-correcting on
    # the new mesh.
    tensor_materials = 'be_bar stress rotation_tensor stretch_tensor stretch_tensor_fbar'
    materials = 'effective_plastic_strain'
  []
[]

[Materials]
  [coalescence]
    type = ADParsedMaterial
    property_name = coal
    material_property_names = 'effective_plastic_strain'
    constant_names = 'c1 c2 c3'
    constant_expressions = '${c1} ${c2} ${c3}'
    expression = ((1-c1)/(1+exp((effective_plastic_strain-c2)/c3))+c1)
    outputs = exodus
    output_properties = 'coal'
  []
  [defgrad]
    type = ComputeDeformationGradient

    output_properties = 'deformation_gradient Fnobar'
    outputs = exodus
  []
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc psic density thermal_conductivity specific_heat'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic} ${rho} ${thermal_conductivity} ${specific_heat}'
  []
  [reg_density]
    type = MaterialADConverter
    ad_props_in = 'density'
    reg_props_out = 'reg_density'
  []
  [nodeg]
    type = NoDegradation
    phase_field = d_corr
    property_name = nodeg
  []
  [g]
    type = PowerDegradationFunction
    phase_field = d_corr
    property_name = g
    expression = (1-d_corr)^p*(1-eta)+eta
    parameter_names = 'p eta'
    parameter_values = '2 1e-6'
  []
  [hencky]
    type = CNHIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d_corr
    degradation_function = g
    decomposition = NONE
  []
  [J2]
    type = LargeDeformationJ2PlasticityBeBar
    phase_field = d_corr
    hardening_model = JC
    relative_tolerance = 1e-08
    output_properties = 'effective_plastic_strain'
    outputs = exodus
  []
  [JC]
    type = JohnsonCookHardening
    T = T
    taylor_quinney_factor = ${Q}
    T0 = 280
    sigma_0 = 1
    reference_plastic_strain = 1
    reference_plastic_strain_rate = 1e-6
    phase_field = d_corr
    degradation_function = g
    A = 791.2
    B = 509.51
    C = 0.014
    n = 0.26
    m = 1.03
    Tm = 1033
    output_properties = 'plastic_heat_generation'
    outputs = exodus
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
    plasticity_model = J2
  []
[]

trans_time = 1.0
final_velocity = 0.2

[Functions]
  [ypull_func]
    type = ParsedFunction
    expression = 'if(t<=trans,v/(2*trans)*t*t,v*t-v*trans/2)'
    symbol_names = 'trans v'
    symbol_values = '${trans_time} ${final_velocity}'
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = 'left top'
    value = 0
    preset = false
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
    preset = false
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = 'front back'
    value = 0
    preset = false
  []
  [ypull]
    type = PresetDisplacement
    variable = disp_y
    boundary = top
    function = ypull_func
    beta = ${newmark_beta}
    velocity = vel_y
    acceleration = accel_y
  []
[]

[Postprocessors]
  [psie_corr_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
  []
  [psip_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psip_active
    use_displaced_mesh = true
  []
  [ep_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = effective_plastic_strain
    use_displaced_mesh = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeIntegrator]
    type = ImplicitEuler
  []

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 50

  [TimeStepper]
    type = ConstantDT
    dt = ${dt}
  []
  [Quadrature]
    order = THIRD
  []

  end_time = ${end_time}

  automatic_scaling = true

  fixed_point_max_its = 25
  fixed_point_rel_tol = 1e-7
  fixed_point_abs_tol = 1e-8
  accept_on_max_fixed_point_iteration = true
  abort_on_solve_fail = true
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = reference_hex8_fbar_out
    use_displaced = false
  []
  [exodusqp]
    type = Exodus
    file_base = reference_hex8_fbar_out_disp
    use_displaced = true
    execute_on = 'FINAL'
  []
  [csv]
    type = CSV
    file_base = reference_hex8_fbar_out
  []
[]
