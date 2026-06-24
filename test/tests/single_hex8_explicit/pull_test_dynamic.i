# Reference run -- single-HEX8 explicit dynamics with central-difference
# time integration (MOOSE's ExplicitMixedOrder).
#
# Companion to pull_test_dynamic_restart.i which exercises the recovery
# pattern (ExplicitMixedOrderRecover seeds the integrator's UDot/UDotDot
# from the dump's nodal vel/accel AuxVariables).
#
# Differences vs the implicit (HHT) variant in single_hex8_incompressible:
#   - No hht_alpha / newmark_beta / newmark_gamma.
#   - No ADInertialForce; CD's integrator handles inertia internally.
#     We add MassMatrix kernels per disp component to assemble the mass
#     matrix into the 'mass' tag.
#   - Standard ADStressDivergenceTensors (no HHT *Recover variant).
#   - No SolutionTensor materials, no recompute_old_stress.
#   - No nl_abs_tol / nl_rel_tol -- CD bypasses the nonlinear solver.
#   - BC is the stock ExplicitFunctionDirichletBC (no spatial variant yet).
#   - vel_*, accel_* are AuxVariables computed via TimeDerivativeAux /
#     SecondTimeDerivativeAux from the integrator's internal time
#     derivatives of disp_*, so they get dumped to exodus for the restart.
#   - Per-QP Jacobian dumped (carries the F-bar mass-matching fix).
#
# Material params: ν = 0.3 (not near-incompressible) to keep the CFL-limited
# dt reasonable; near-incompressible cases need much smaller dt.

E = 1.0
nu = 0.3
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 1.0

# CFL-limited dt.  Wave speed c = sqrt((K + 4G/3)/ρ) ≈ 1.16 for ν=0.3.
# Element size 1 → dt_CFL ≤ ~0.86.  Use dt well below that.
dt = 0.001
end_time = 1
pull_amount = 1

out_dir = outputs
tag = ''
dump_time = 0.5
end_time_for_ramp = ${end_time}

xfix_bnd = 'left'
yfix_bnd = 'bottom'
zfix_bnd = 'front back'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
  # HEX8_3rd for QP indexing consistent with the recovery infrastructure.
  element = HEX8_3rd
[]

[Problem]
  type = FEProblem
  # Mass matrix lives in this tag; ExplicitMixedOrder reads from it.
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
    elem_type = HEX8
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
  # Phase-field stub, kept zero; the recovery infrastructure expects it.
  [d_corr]
  []
  # Original undeformed nodal x-coordinate, frozen at t=0.  Dumped so the
  # restart can recover it; not used by the simple uniform BC here.
  [X0]
    family = LAGRANGE
    order = FIRST
    [InitialCondition]
      type = FunctionIC
      function = 'x'
    []
  []

  # ----- CD's internal vel / accel mirrored to AuxVariables for the dump --
  # The restart's ExplicitMixedOrderRecover reads these from the dump's
  # exodus and writes them directly into its integrator state (UDot/UDotDot).
  [vel_x]
  []
  [vel_y]
  []
  [vel_z]
  []
  [accel_x]
  []
  [accel_y]
  []
  [accel_z]
  []
[]

[AuxKernels]
  # Nodal copy of the CD integrator's internal velocity / acceleration into
  # AuxVariables so they can be dumped to exodus.  CoupledTimeDerivativeAux
  # (raccoon-side) uses coupledDot / coupledDotDot which work for both
  # nodal and elemental targets, unlike MOOSE's TimeDerivativeAux which is
  # elemental-only.
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
  # Internal force (no HHT, no inertia term in the kernel -- the integrator
  # handles inertia separately via the MassMatrix kernels + central-diff
  # solve).
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

  # ----- MassMatrix kernels (one per disp component) -----------------------
  # ExplicitMixedOrder assembles these into the 'mass' matrix tag, then
  # row-sum-lumps and inverts the diagonal.  density = constant rho here
  # (no recovery-side mass adjustment on the reference).
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
    # Linear ramp.  No spatial dependence (top face moves uniformly).
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
    # Explicit-integrator-aware Dirichlet BC.  Reads the lumped mass diag
    # vector from the integrator and computes a residual that drives disp_y
    # at the top face to pull_func(t) after the CD update.  Spatial variant
    # (coupled_x = X0) deferred -- this uses node->point() directly.
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
    # K, G stay AD because the elasticity model is AD.
    type = ADGenericConstantMaterial
    prop_names = 'K G'
    prop_values = '${K} ${G}'
  []
  [density]
    # density must be non-AD: MassMatrix declares a non-AD density getter.
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
  [disp_y_top_int]
    # Integrated top-face disp_y -- proxy for how far the loading has pushed.
    type = SideAverageValue
    variable = disp_y
    boundary = top
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

  end_time = ${end_time}
  # No nl_*_tol / nl_max_its -- ExplicitMixedOrder bypasses the nonlinear solver.
[]

[RecoverVariables]
  [rec]
    # Approach C requires raw stretch_tensor (U from polar decomp of F_raw).
    # libmesh's exodus reader shadows the shorter name when both `stretch_tensor`
    # and `stretch_tensor_fbar` are present, so omit fbar to make sure raw
    # stretch_tensor is what ends up in the file.
    tensor_materials = 'stress rotation_tensor stretch_tensor stretch_tensor_fbar'
    # Per-QP Jacobian: used by both CDG approach C (denom_C formula in the
    # restart's CDG) and the restart's adj_density_rec.
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
    # Recovery dump: one slice at dump_time, displaced mesh, so the restart
    # can FileMeshGenerator-load it as its undisplaced reference.
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
