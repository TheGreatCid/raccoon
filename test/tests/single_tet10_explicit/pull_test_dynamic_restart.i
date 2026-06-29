# Restart leg of the single-TET10 CD recovery test.  Pairs with
# pull_test_dynamic.i (reference).
#
# Mirrors single_hex8_explicit/pull_test_dynamic_restart.i with the
# TET10-specific adaptations:
#   - Variables/AuxVariables order = SECOND family = LAGRANGE so
#     SolutionIC (for X0/d_corr) AND the postprocessor reads see the
#     mid-edge node values.
#   - element = TET10_4th in GlobalParams and in SolutionReal*.
#   - Quadrature order = FOURTH.

E = 1.0
nu = 0.3
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 1.0

dt = 0.0001
start_time = 0.025
end_time = 0.05
pull_amount = 0.005
end_time_for_ramp = 0.05
# Pulse parameters must match the reference's exactly.  See pull_test_dynamic.i
# for rationale.
pulse_amplitude = 0.005
pulse_center = 0.01
pulse_width = 0.02

out_dir = outputs
tag = ''
output_tag = ${tag}
recover_file = ${out_dir}/pull_test_out_disp${tag}.e

xfix_bnd = 'left'
yfix_bnd = 'bottom'
zfix_bnd = 'front back'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
  element = TET10_4th
[]

[Problem]
  type = FEProblem
  extra_tag_matrices = 'mass'
[]

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = ${recover_file}
  []
[]

[UserObjects]
  [epsol]
    type = SolutionUserObjectQP
    mesh = ${recover_file}
    # Approach C requires raw stretch_tensor.
    tensor_materials = 'stress stretch_tensor rotation_tensor'
    materials = 'Jacobian'
    system_variables = 'd_corr X0 accel_x accel_y accel_z vel_x vel_y vel_z'
    use_displaced_mesh = true
    execute_on = 'INITIAL'
    timestep = LATEST
  []
[]

[Variables]
  [disp_x]
    order = SECOND
  []
  [disp_y]
    order = SECOND
  []
  [disp_z]
    order = SECOND
  []
[]

[AuxVariables]
  [d_corr]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = d_corr
      variable = d_corr
      solution_uo = epsol
    []
  []
  [X0]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = X0
      variable = X0
      solution_uo = epsol
    []
  []

  # vel_*, accel_* on restart: order = SECOND family = LAGRANGE.  The
  # integrator state itself is seeded by ExplicitMixedOrderRecover via
  # direct vector writes; these AuxVars are repopulated each step by
  # CoupledTimeDerivativeAux for diagnostics / downstream output.
  [vel_x]
    order = SECOND
    family = LAGRANGE
  []
  [vel_y]
    order = SECOND
    family = LAGRANGE
  []
  [vel_z]
    order = SECOND
    family = LAGRANGE
  []
  [accel_x]
    order = SECOND
    family = LAGRANGE
  []
  [accel_y]
    order = SECOND
    family = LAGRANGE
  []
  [accel_z]
    order = SECOND
    family = LAGRANGE
  []
[]

[AuxKernels]
  [vel_x_aux]
    type = CoupledTimeDerivativeAux
    variable = vel_x
    v = disp_x
    execute_on = 'timestep_end'
  []
  [vel_y_aux]
    type = CoupledTimeDerivativeAux
    variable = vel_y
    v = disp_y
    execute_on = 'timestep_end'
  []
  [vel_z_aux]
    type = CoupledTimeDerivativeAux
    variable = vel_z
    v = disp_z
    execute_on = 'timestep_end'
  []
  [accel_x_aux]
    type = CoupledTimeDerivativeAux
    variable = accel_x
    v = disp_x
    order_derivative = SECOND
    execute_on = 'timestep_end'
  []
  [accel_y_aux]
    type = CoupledTimeDerivativeAux
    variable = accel_y
    v = disp_y
    order_derivative = SECOND
    execute_on = 'timestep_end'
  []
  [accel_z_aux]
    type = CoupledTimeDerivativeAux
    variable = accel_z
    v = disp_z
    order_derivative = SECOND
    execute_on = 'timestep_end'
  []
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

  [mass_x]
    type = MassMatrix
    variable = disp_x
    density = adj_density
    matrix_tags = 'mass'
  []
  [mass_y]
    type = MassMatrix
    variable = disp_y
    density = adj_density
    matrix_tags = 'mass'
  []
  [mass_z]
    type = MassMatrix
    variable = disp_z
    density = adj_density
    matrix_tags = 'mass'
  []
[]

[Functions]
  [ypull_func_restart]
    # INCREMENTAL top-y displacement past start_time.  Both the ramp AND the
    # pulse are subtracted at start_time so the increment is 0 at the
    # restart's first step (matching the dump's stored disp_y at the top).
    type = ParsedFunction
    expression = 'pull_amount * ((t - start_time) / end_time_for_ramp) + pulse_amplitude * (exp(-((t-pulse_center)*(t-pulse_center))/(pulse_width*pulse_width)) - exp(-((start_time-pulse_center)*(start_time-pulse_center))/(pulse_width*pulse_width)))'
    symbol_names = 'pull_amount start_time end_time_for_ramp pulse_amplitude pulse_center pulse_width'
    symbol_values = '${pull_amount} ${start_time} ${end_time_for_ramp} ${pulse_amplitude} ${pulse_center} ${pulse_width}'
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = '${xfix_bnd}'
    value = 0
    preset = true
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = '${yfix_bnd}'
    value = 0
    preset = true
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = '${zfix_bnd}'
    value = 0
    preset = true
  []
  [ypull]
    type = ExplicitFunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ypull_func_restart
  []
[]

[Materials]
  [defgrad]
    type = ComputeDeformationGradient
    recover = true
    solution = epsol
    recover_mode = polar_decomposition
    volumetric_locking_correction = true
    recover_apply_fbar_to_total = true
    output_properties = 'deformation_gradient Fnobar'
    outputs = exodus
  []
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [rho_const]
    type = GenericConstantMaterial
    prop_names = 'rho_const'
    prop_values = '${rho}'
  []
  # Geometric-J adjusted density (det of _F_store_noFbar, which after polar-
  # decomp recovery equals the dump-deformed mesh's local Jacobian).  Same
  # rationale as in single_hex8_explicit.
  [adj_density]
    type = StrainAdjustedDensityCustom
    strain_free_density = rho_const
    base_name = 'adj'
  []
  [nodeg]
    type = NoDegradation
    phase_field = d_corr
    property_name = nodeg
  []
  [hencky]
    type = CNHIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d_corr
    degradation_function = nodeg
    decomposition = NONE
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
  []
[]

[Postprocessors]
  [psie_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
  []
  [psie_at_init]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
    execute_on = 'INITIAL'
  []
  [vel_y_max]
    type = NodalExtremeValue
    variable = vel_y
    value_type = max_abs
  []
  [accel_y_max]
    type = NodalExtremeValue
    variable = accel_y
    value_type = max_abs
  []
[]

[Executioner]
  type = Transient

  [TimeIntegrator]
    type = ExplicitMixedOrderRecover
    mass_matrix_tag = 'mass'
    second_order_vars = 'disp_x disp_y disp_z'
    solution   = epsol
    disp_vars  = 'disp_x disp_y disp_z'
    vel_vars   = 'vel_x vel_y vel_z'
    accel_vars = 'accel_x accel_y accel_z'
    use_constant_mass = true
  []

  [TimeStepper]
    type = ConstantDT
    dt = ${dt}
  []

  [Quadrature]
    order = FOURTH
  []

  start_time = ${start_time}
  end_time = ${end_time}
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = ${out_dir}/pull_test_restart_out${output_tag}
    use_displaced = false
  []
  [csv]
    type = CSV
    file_base = ${out_dir}/pull_test_restart_out${output_tag}
  []
[]
