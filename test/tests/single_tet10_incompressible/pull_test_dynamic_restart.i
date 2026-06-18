# Restart leg of the single-TET10 incompressible bending recovery test
# (DYNAMIC variant).
# Pairs with pull_test_dynamic.i (reference).
#
# Same recovery / X0 / MatchedValueBC structure as pull_test_restart.i but
# with Newmark/HHT inertia added.  accel_*, vel_* are recovered from the
# dump via SolutionIC so the body's momentum picks up where the reference
# left it at start_time.  NewmarkAccelAux / NewmarkVelAux keep the state
# evolving from there.
#
# Same physics suite as pull_test.i (elastic + quasi-static, near-incompressible
# CNH bending), but the mesh comes from FileMeshGenerator(file = recover_file)
# so the body starts in the reference's deformed configuration at dump_time.
# The kinematic state (R, U or U_bar) is recovered through a SolutionUserObjectQP
# and seeded into ComputeDeformationGradient(recover = true).
#
# Two recovery methods are supported via the CLI:
#   Approach A (NEW, default):   recover U_bar directly,
#                                Materials/defgrad/recover_apply_fbar_to_U=false
#   Approach B (OLD):            recover raw U and apply F-bar averaging at
#                                INITIAL via recover_apply_fbar_to_U=true.

E = 1.0
nu = 0.4999
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 1.0

# HHT-alpha / Newmark-beta time-integration parameters.
hht_alpha = -0.25
newmark_beta = '${fparse (1-hht_alpha)^2/4}'
newmark_gamma = '${fparse 1/2-hht_alpha}'

# Time stepping: must match the reference.i values.
dt = 0.0125
start_time = 0.5
end_time = 1.0
pull_amount = 0.2
# The reference's end_time at which the loading would saturate -- used to
# scale ypull_func / ypull_func_restart so the BC ramp is the same in both
# legs.  Always the original reference's end_time (= 1.0 by default).
end_time_for_ramp = 1.0

# Pulse + ramp_time (must match the reference's values).  Pulse splits into
# bending (linear in X0) + uniform compression.
pulse_amplitude_bend = 0.01
pulse_amplitude_compress = 0.01
pulse_center = 0.05
pulse_width = 0.025
ramp_time = 0.1

# Output / recovery wiring.
out_dir = outputs
tag = ''
output_tag = ${tag}
recover_file = ${out_dir}/pull_test_out_disp${tag}.e
# LATEST is correct because the reference dumps exactly one slice at dump_time
# via sync_times+sync_only.
recover_timestep = LATEST

# BC face lists.  Must match the reference run for the recovered state to
# remain in equilibrium.
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
    # nu_sweep_study.sh overrides this list per method so the restart asks for
    # only the stretch tensor variant the reference dumped (libmesh's exodus
    # reader shadows the shorter name when both prefixes are present).
    tensor_materials = 'stress stretch_tensor_fbar rotation_tensor'
    # Jacobian: per-QP J_raw for adj_density_rec.  See
    # single_hex8_incompressible/dynamic_recovery_fix.pdf.
    materials = 'Jacobian'
    # Recover X0 (frozen undeformed nodal x-coords) so the BC can use the
    # ORIGINAL x rather than the deformed-mesh coordx, and d_corr (phase
    # field stub).
    system_variables = 'd_corr X0 accel_x accel_y accel_z vel_x vel_y vel_z'
    nodal_variable_order = SECOND
    use_displaced_mesh = true
    execute_on = 'INITIAL'
    timestep = ${recover_timestep}
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
    # Phase-field stub.  Recovered (will be ~0 if the reference left it at 0).
    [InitialCondition]
      type = SolutionIC
      from_variable = d_corr
      variable = d_corr
      solution_uo = epsol
    []
  []
  [X0]
    # Recover the reference's frozen undeformed nodal x-coords.  Without
    # this, the restart's BC would use the deformed mesh's coordx for x in
    # the bending function, which has already laterally contracted from the
    # reference's t=dump_time state.
    family = LAGRANGE
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = X0
      variable = X0
      solution_uo = epsol
    []
  []
  [target_uy]
    # BC target for disp_y on the top, computed from X0 instead of the
    # deformed mesh's x.  ParsedAux below assembles
    #   target_uy = pull_amount * X0 * (t - start_time) / end_time_for_ramp
    # i.e. the INCREMENTAL displacement past start_time, which is what the
    # restart needs to add on top of the recovered deformed configuration.
    family = LAGRANGE
    order = SECOND
  []

  # ----- Newmark / HHT dynamic state (recovered from dump) ------------------
  # CRITICAL: order = SECOND must match the reference's accel_*/vel_*
  # AuxVariable order.  Without this, SolutionIC reads from a SECOND-order
  # source into a default FIRST-order target, dropping mid-edge node values
  # and producing wildly wrong accel/vel at the recovery boundary.
  [accel_x]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_x
      variable = accel_x
      solution_uo = epsol
    []
  []
  [vel_x]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_x
      variable = vel_x
      solution_uo = epsol
    []
  []
  [accel_y]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_y
      variable = accel_y
      solution_uo = epsol
    []
  []
  [vel_y]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_y
      variable = vel_y
      solution_uo = epsol
    []
  []
  [accel_z]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_z
      variable = accel_z
      solution_uo = epsol
    []
  []
  [vel_z]
    order = SECOND
    family = LAGRANGE
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_z
      variable = vel_z
      solution_uo = epsol
    []
  []
[]

[AuxKernels]
  [compute_target_uy]
    # Diagnostic AuxVariable (BC now consumes ypull_func_restart directly).
    type = ParsedAux
    variable = target_uy
    coupled_variables = 'X0'
    constant_names       = 'pull_amount start_time end_time_for_ramp pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    constant_expressions = '${pull_amount} ${start_time} ${end_time_for_ramp} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'
    expression = 'pull_amount * X0 * (t - start_time) / end_time_for_ramp + (pulse_amplitude_bend * X0 + pulse_amplitude_compress) * (exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width)) - exp(-(start_time-pulse_center)*(start_time-pulse_center)/(pulse_width*pulse_width)))'
    use_xyzt = true
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR'
  []

  # Newmark accel/vel updates.
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
  # Recover-aware HHT stress divergence: at the first restart step the regular
  # MOOSE _stress_old is zero (no prior state on the restart side), which drops
  # the -alpha * sigma_old * grad_test term from the HHT-alpha residual and
  # produces a 25% missing force on step 1.  The *Recover variant reads
  # sigma_old from the SolutionUserObject (= recovered reference stress at
  # dump_time) for t_step=1, then transitions to MOOSE-tracked _stress_old at
  # t_step >= 2.  Requires the two SolutionTensor materials below.
  [x]
    type = ADDynamicStressDivergenceTensorsRecover
    variable = disp_x
    component = 0
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    solution = epsol
  []
  [y]
    type = ADDynamicStressDivergenceTensorsRecover
    variable = disp_y
    component = 1
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    solution = epsol
  []
  [z]
    type = ADDynamicStressDivergenceTensorsRecover
    variable = disp_z
    component = 2
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    solution = epsol
  []

  # ----- Inertia (recovered velocity/acceleration carry through) ------------
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    density = adj_density_rec
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_x
    acceleration = accel_x
  []
  [inertia_y]
    type = ADInertialForce
    variable = disp_y
    density = adj_density_rec
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_y
    acceleration = accel_y
  []
  [inertia_z]
    type = ADInertialForce
    variable = disp_z
    density = adj_density_rec
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_z
    acceleration = accel_z
  []
[]

[Functions]
  [ypull_func]
    # Absolute reference top-y displacement: pull_amount * (t/end_time_for_ramp) * x
    # -- linear-in-x bending as in pull_test.i.
    type = ParsedFunction
    expression = 'pull_amount * (t/end_time_for_ramp) * x'
    symbol_names = 'pull_amount end_time_for_ramp'
    symbol_values = '${pull_amount} ${end_time_for_ramp}'
  []
  [ypull_func_restart]
    # INCREMENTAL top-y displacement on the restart's deformed-config mesh.
    # Assumes start_time > ramp_time so the ramp_time/2 offset cancels in
    # the incremental form.  Pulse term subtracts the pulse value at
    # start_time so the BC stays continuous across the recovery boundary.
    type = ParsedFunction
    expression = 'pull_amount * ((t - start_time)/end_time_for_ramp) * x + (pulse_amplitude_bend * x + pulse_amplitude_compress) * (exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width)) - exp(-(start_time-pulse_center)*(start_time-pulse_center)/(pulse_width*pulse_width)))'
    symbol_names = 'pull_amount start_time end_time_for_ramp pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    symbol_values = '${pull_amount} ${start_time} ${end_time_for_ramp} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = '${xfix_bnd}'
    value = 0
    preset = false
  []
  [yfix]
    type = DirichletBC
    variable = disp_y
    boundary = '${yfix_bnd}'
    value = 0
    preset = false
  []
  [zfix]
    type = DirichletBC
    variable = disp_z
    boundary = '${zfix_bnd}'
    value = 0
    preset = false
  []
  [ypull]
    # PresetDisplacementSpatial: raccoon's spatially-aware variant of
    # MOOSE's PresetDisplacement.  Routes the function's spatial argument
    # through coupled_x = X0 (the recovered AuxVariable carrying the original
    # undeformed x of each node), so the spatial gradient comes from the
    # original undeformed mesh, not from the dump-deformed mesh's coordx.
    type = PresetDisplacementSpatial
    variable = disp_y
    boundary = top
    function = ypull_func_restart
    coupled_x = X0
    beta = ${newmark_beta}
    velocity = vel_y
    acceleration = accel_y
  []
[]

[Materials]
  # ----- HHT recover-side stress history --------------------------------------
  # ADDynamicStressDivergenceTensorsRecover requires `stress_sol` (the recovered
  # reference stress at dump_time, seeded via the SolutionUserObject) and a
  # `stress_old_store_sol` placeholder (only used multiplied by zeta, which is
  # 0 here, so the identity default is harmless).
  [stress_sol]
    type = SolutionTensor
    solution = epsol
    tensor_name = stress
  []
  [stress_old_sol]
    type = ADGenericConstantRankTwoTensor
    tensor_name = stress_old_store_sol
    tensor_values = '1 0 0 0 1 0 0 0 1'
  []
  [defgrad]
    type = ComputeDeformationGradient
    recover = true
    solution = epsol
    recover_mode = polar_decomposition
    # nu_sweep_study.sh overrides recover_apply_fbar_to_U to toggle approach
    # A (default false, consume U_bar) vs approach B (true, rebuild F_bar
    # from raw U via averaging at INITIAL).
    output_properties = 'deformation_gradient Fnobar'
    outputs = exodus
  []
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G density'
    prop_values = '${K} ${G} ${rho}'
  []
  [dens]
    # Kept for diagnostics (det(_F_NoFbar) which in approach A is element-
    # constant J_bar).  Inertia kernels consume adj_density_rec below, which
    # uses the per-QP J_raw recovered from the reference's dump.  See
    # single_hex8_incompressible/dynamic_recovery_fix.pdf.
    type = ADStrainAdjustedDensityCustom
    strain_free_density = density
    base_name = 'adj'
  []
  # ---- Per-QP J_raw recovery for exact inertia mass matching ----
  [recovered_J]
    type = SolutionReal
    solution = epsol
    mat_name = Jacobian
    element = TET10_4th
  []
  [adj_density_rec]
    type = ADParsedMaterial
    property_name = adj_density_rec
    material_property_names = 'Jacobian_sol'
    constant_names = 'rho'
    constant_expressions = '${rho}'
    expression = 'rho / Jacobian_sol'
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

  # F-bar diagnostic (same definitions as pull_test.i).
  [J_F]
    type = ADRankTwoInvariant
    rank_two_tensor = Fnobar
    invariant = ThirdInvariant
    property_name = J_F
    outputs = exodus
  []
  [J_Fbar]
    type = ADRankTwoInvariant
    rank_two_tensor = deformation_gradient
    invariant = ThirdInvariant
    property_name = J_Fbar
    outputs = exodus
  []
  [fbar_correction]
    type = ADParsedMaterial
    property_name = fbar_correction
    material_property_names = 'J_F J_Fbar'
    expression = abs(J_F/J_Fbar-1)
    outputs = exodus
  []
  [fbar_pressure_correction]
    type = ADParsedMaterial
    property_name = fbar_pressure_correction
    material_property_names = 'J_F J_Fbar K'
    expression = K*abs(J_Fbar-J_F)
    outputs = exodus
  []
[]

[Postprocessors]
  [psie_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
  []
  # Diagnostics: state at INITIAL of restart (should match ref's state at t=dump_time)
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
  [J_F_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = J_F
    use_displaced_mesh = false
  []
  [J_Fbar_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = J_Fbar
    use_displaced_mesh = false
  []
  [fbar_correction_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = fbar_correction
    use_displaced_mesh = false
  []
  [fbar_correction_max]
    type = ADElementExtremeMaterialProperty
    mat_prop = fbar_correction
    value_type = max
  []
  [fbar_pressure_correction_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = fbar_pressure_correction
    use_displaced_mesh = false
  []
  [fbar_pressure_correction_max]
    type = ADElementExtremeMaterialProperty
    mat_prop = fbar_pressure_correction
    value_type = max
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  [TimeIntegrator]
    type = ImplicitEuler
  []

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 50

  [TimeStepper]
    type = ConstantDT
    dt = ${dt}
  []
  [Quadrature]
    order = FOURTH
  []

  start_time = ${start_time}
  end_time = ${end_time}
  automatic_scaling = true
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
