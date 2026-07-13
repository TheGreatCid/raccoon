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
dt =0.5 #0.5
end_time = 400
pull_amount = 0.9


# BC ramp denominator and Gaussian pulse params (see pull_test_dynamic.i in
# single_hex8_incompressible for the rationale).
end_time_for_ramp = 200
pulse_amplitude_bend = 0
pulse_amplitude_compress = 0.8#0.4
# Quadratic-then-linear startup window so f''(t) is continuous at t=0,
# required by PresetDisplacementSpatial.
ramp_time = 1
pulse_center = 50
pulse_width = 10
xfix_bnd = 'left_ns'
yfix_bnd = 'bottom_ns'
zfix_bnd = 'front_ns back_ns'

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
    type = FileMesh
    construct_node_list_from_side_list=False
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
    [frob_time]
    order=CONSTANT
    family=MONOMIAL
  []
  [psie_old]
    order=CONSTANT
    family=MONOMIAL
  []
  [Y0]
  []
  [d]
  []
  [sizing]

  []
 [qual_frob]
    order = CONSTANT
    family = MONOMIAL
  []

  # Phase-field stub, kept zero; the recovery infrastructure expects it.
  [d_corr]
  []
  # Original undeformed nodal x-coordinate, frozen at t=0.  Dumped so the
  # restart can recover it; not used by the simple uniform BC here.
  [X0]
    family = LAGRANGE
    order = FIRST
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
  [kinetic_energy]
    order=CONSTANT
    family=MONOMIAL
  []
  [kinetic_energy_lf]
    order=CONSTANT
    family=MONOMIAL
  []
[]

[AuxKernels]
  [psie_old]
    type = ADMaterialRealAux
    property=psie_active
    variable = psie_old 
    execute_on = 'TIMESTEP_BEGIN'
  []
    [qual_frob]
    type = MaterialRealAux
    property = Frobenius_norm
    variable = qual_frob
    execute_on = 'TIMESTEP_END'
  []
  [ke_lf]
    type = LeapfrogKineticEnergyAux
    velocity_x = vel_x
    velocity_y = vel_y
    velocity_z = vel_z
    density = density
    variable = kinetic_energy_lf
    execute_on  = 'timestep_end'
  []
    [frob_time]
    type = ParsedAux
    variable = frob_time
    expression = 'if(t<360,0,qual_frob)'
    use_xyzt=true
    coupled_variables = qual_frob
  []

  [Ke]
    type = KineticEnergyAux
    newmark_velocity_x = vel_x
    newmark_velocity_y = vel_y
    newmark_velocity_z = vel_z
    density = density
    variable = kinetic_energy
  []

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
 ## [vel_x_aux]
 #   type = TestNewmarkTI
 #   variable = vel_x
 #   displacement = disp_x
 #   execute_on = 'linear timestep_begin timestep_end'
 # []
 # [vel_y_aux]
 #   type = TestNewmarkTI
 #   variable = vel_y
 #   displacement = disp_y
 #   execute_on = 'linear timestep_begin timestep_end'
 # []
 # [vel_z_aux]
 #   type = TestNewmarkTI
 #   variable = vel_z
 #   displacement = disp_z
 #   execute_on = 'linear timestep_begin timestep_end'
 # []
 # [accel_x_aux]
 #   type = TestNewmarkTI
 #   variable = accel_x
 #   displacement = disp_x
 #   first='false'
 #   execute_on = 'linear timestep_begin timestep_end'
 # []
 # [accel_y_aux]
 #   type = TestNewmarkTI
 #   variable = accel_y
 #   displacement = disp_y
 #   first='false'
 #   execute_on = 'linear timestep_begin timestep_end'
 # []
 # [accel_z_aux]
 #   type = TestNewmarkTI
 #   variable = accel_z
 #   displacement = disp_z
 #   first='false'
 #   execute_on = 'linear timestep_begin timestep_end'
 # []
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
    [cent]
    type = ParsedFunction
    expression = '1'
    #symbol_names = 'minval maxval thick mindist dist'
    #symbol_values = '${minval} ${maxval} ${thick} ${mindist} dist'
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
    [frob]
    type = ElementExtremeValue
    variable = frob_time
    value_type = max
    execute_on = 'TIMESTEP_END'
  []
  [kinetic_energy]

    type = ElementIntegralVariablePostprocessor
    variable = kinetic_energy
    use_displaced_mesh = true
  []
   [kinetic_energy_int]
    type = ElementIntegralVariablePostprocessor
    variable = kinetic_energy_lf
    use_displaced_mesh = true
  []
  [psie_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true

    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [psie_old_int]
    type = ElementIntegralVariablePostprocessor
    variable = psie_old
    use_displaced_mesh = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [total_energy]
    type = ParsedPostprocessor
    expression = 'psie_active_int + kinetic_energy_int'
    pp_names = 'psie_active_int kinetic_energy_int'
    #use_displaced_mesh=true
    execute_on = 'TIMESTEP_END'
  []
  # KE in the integrator's OWN lumped-mass norm, 1/2 sum_i m_i (u_dot_i)^2.
  # Comparing this between the reference and recover legs isolates the lumped-
  # mass reconstruction error (recovered adj_density/det(F) vs the reference's
  # rho0-on-Omega0 mass), since the recovered velocity matches to ~1e-4.
  [ke_lumped]
    type = LumpedKineticEnergy
    execute_on = 'TIMESTEP_END'
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

[UserObjects]
  [Terminator]
    type = Terminator
  #expression = 'frob > 0.0001'
  expression = '1<0'
    fail_mode = HARD
    error_level = ERROR
    message = 'MESH FLAG'
    force_postaux = false
    execute_on = MULTIAPP_FIXED_POINT_BEGIN
  []
[]

[Outputs]

  print_linear_residuals = false
time_step_interval = 1
    [exodusqp]
    type = Exodus
    use_displaced = false
#    execute_on = 'FINAL'
  []
  [exodus]
    type = Exodus
  []
  [csv]
    type = CSV
    file_base = ${csv_path}
    time_step_interval=1
  []
[]
