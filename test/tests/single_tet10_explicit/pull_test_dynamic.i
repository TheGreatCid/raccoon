# Reference run -- single-TET10 (one HEX8 split into ~6 TET10 elements) with
# central-difference explicit dynamics.  Companion to
# pull_test_dynamic_restart.i.
#
# Differences vs single_hex8_explicit/pull_test_dynamic.i:
#   - Mesh element type: TET10 (HEX8 split).
#   - GlobalParams element = TET10_4th (4 QPs/tet; raccoon's QP-index
#     mapping for the recovery infrastructure).
#   - Variables disp_x/y/z declared order = SECOND (TET10 carries DOFs
#     at corner + mid-edge nodes).
#   - AuxVariables (vel_*, accel_*, X0, d_corr) all explicitly
#     order = SECOND family = LAGRANGE.  This matters for the dump-to-
#     restart round-trip: if the restart-side AuxVars default to FIRST,
#     SolutionIC silently drops mid-edge node values and breaks recovery
#     (this is the same gotcha documented for the implicit TET10 case).
#   - Quadrature order = FOURTH (matches TET10_4th).
#   - dt smaller for CFL stability on the smaller tet elements.

E = 1.0
nu = 0.3
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 1.0

# CFL-limited dt.  TET10 elements after the HEX8 split are smaller than
# the parent hex; the smallest sub-tet edge is roughly the cube edge / 2
# (in the worst case for a 5- or 6-tet split).  Wave speed
# c = sqrt((K + 4G/3)/rho) ~ 1.16 for nu=0.3.  Element size ~ 0.5 ->
# dt_CFL ~ 0.43.  Use dt = 1e-4 to be safely under (TET10's higher-order
# basis can have more restrictive CFL than the geometric bound suggests).
dt = 0.0001
# Short test window -- end_time = 0.05 keeps the run under ~30s and the
# trajectory divergence below ~1% over the post-dump window.  CD on TET10
# accumulates per-step error from the same lumped-mass-mismatch source as
# the HEX8 test (J_recovered vs J_mesh_geometric not cancelling exactly in
# (rho/J) * _JxW), but the per-step rate is larger on TET10 because of the
# higher mode density of the 24-element mesh.  Longer runs amplify this
# drift; for recovery-infrastructure validation, the short window is
# sufficient.
end_time = 0.05
pull_amount = 0.1

out_dir = outputs
tag = ''
dump_time = 0.025
end_time_for_ramp = ${end_time}

xfix_bnd = 'left'
yfix_bnd = 'bottom'
zfix_bnd = 'front back'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
  # TET10_4th matches FOURTH-order quadrature (4 QPs / tet).  Required by
  # the recovery infrastructure for consistent QP indexing.
  element = TET10_4th
[]

[Problem]
  type = FEProblem
  extra_tag_matrices = 'mass'
[]

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1
    elem_type = TET10
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
  []
  [X0]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = FunctionIC
      function = 'x'
    []
  []

  # vel_*, accel_* mirror CD's internal state for the dump.  ORDER = SECOND
  # FAMILY = LAGRANGE so mid-edge node values get stored to exodus; the
  # restart's directValue(node, var_name) lookup then sees the right values
  # at every node.
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
    density = density
    matrix_tags = 'mass'
  []
  [mass_y]
    type = MassMatrix
    variable = disp_y
    density = density
    matrix_tags = 'mass'
  []
  [mass_z]
    type = MassMatrix
    variable = disp_z
    density = density
    matrix_tags = 'mass'
  []
[]

[Functions]
  [pull_func]
    type = ParsedFunction
    expression = 'pull_amount * (t / end_time_for_ramp)'
    symbol_names = 'pull_amount end_time_for_ramp'
    symbol_values = '${pull_amount} ${end_time_for_ramp}'
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
    function = pull_func
  []
[]

[Materials]
  [defgrad]
    type = ComputeDeformationGradient
    output_properties = 'deformation_gradient Fnobar'
    outputs = exodus
  []
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '${rho}'
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
    type = ExplicitMixedOrder
    mass_matrix_tag = 'mass'
    second_order_vars = 'disp_x disp_y disp_z'
    use_constant_mass = true
  []

  [TimeStepper]
    type = ConstantDT
    dt = ${dt}
  []

  [Quadrature]
    order = FOURTH
  []

  end_time = ${end_time}
[]

[RecoverVariables]
  [rec]
    # Approach C: dump raw stretch_tensor (not stretch_tensor_fbar) so the
    # restart's CDG can reconstruct F_dump_raw and apply F-bar to the total
    # F at every step.
    tensor_materials = 'stress rotation_tensor stretch_tensor'
    materials = 'Jacobian'
  []
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = ${out_dir}/pull_test_out${tag}
    use_displaced = false
  []
  [exodusqp]
    type = Exodus
    file_base = ${out_dir}/pull_test_out_disp${tag}
    use_displaced = true
    sync_times = '${dump_time}'
    sync_only = true
  []
  [csv]
    type = CSV
    file_base = ${out_dir}/pull_test_out${tag}
  []
[]
