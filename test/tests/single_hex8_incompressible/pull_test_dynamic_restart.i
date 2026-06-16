# Restart leg of the single-HEX8 incompressible bending recovery test
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
hht_alpha = -0.1 #-0.25
newmark_beta = '${fparse (1-hht_alpha)^2/4}'
newmark_gamma = '${fparse 1/2-hht_alpha}'

# Time stepping: must match the reference.i values.
dt = 0.0125
start_time = 0.5
end_time = 1.0
pull_amount = 0.2

# Gaussian pulse on the top-face BC (must match pull_test_dynamic.i for the
# reference and restart trajectories to coincide).  The restart's incremental
# BC subtracts the pulse value at start_time so the BC stays continuous
# across the recovery boundary even if the dump happens during the pulse.
# Pulse split (must match the reference's split -- see pull_test_dynamic.i for
# rationale).
pulse_amplitude_bend = 0.01
pulse_amplitude_compress = 0.01
pulse_center = 0.05
pulse_width = 0.025
# The reference's end_time at which the loading would saturate -- used to
# scale ypull_func / ypull_func_restart so the BC ramp is the same in both
# legs.  Always the original reference's end_time (= 1.0 by default).
end_time_for_ramp = 1.0

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
  element = HEX8_3rd
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
    materials = ''
    # Recover X0 (frozen undeformed nodal x-coords) so the BC can use the
    # ORIGINAL x rather than the deformed-mesh coordx, and d_corr (phase
    # field stub).
    system_variables = 'd_corr X0 accel_x accel_y accel_z vel_x vel_y vel_z'
    use_displaced_mesh = true
    execute_on = 'INITIAL'
    timestep = ${recover_timestep}
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
  [d_corr]
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
    order = FIRST
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
    order = FIRST
  []

  # ----- Newmark / HHT dynamic state (recovered from dump) ------------------
  [accel_x]
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_x
      variable = accel_x
      solution_uo = epsol
    []
  []
  [vel_x]
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_x
      variable = vel_x
      solution_uo = epsol
    []
  []
  [accel_y]
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_y
      variable = accel_y
      solution_uo = epsol
    []
  []
  [vel_y]
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_y
      variable = vel_y
      solution_uo = epsol
    []
  []
  [accel_z]
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_z
      variable = accel_z
      solution_uo = epsol
    []
  []
  [vel_z]
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
    type = ParsedAux
    variable = target_uy
    coupled_variables = 'X0'
    constant_names = 'pull_amount start_time end_time_for_ramp pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
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
    density = adj_density
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
    density = adj_density
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
    density = adj_density
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
  [X0_sol]
    type = SolutionFunction
    from_variable = X0
    solution = epsol
  []
  [ypull_func_restart]
    # INCREMENTAL top-y displacement (= reference's pull_func(t) minus its
    # value at start_time) for the post-ramp era.  Consumed by the
    # PresetDisplacementSpatial BC below; the BC reroutes the function's `x`
    # argument to the recovered X0 AuxVariable via coupled_x, so the spatial
    # gradient is the ORIGINAL undeformed x and not the dump-deformed mesh's x.
    #
    # Assumes start_time > ramp_time of the reference run (typical for the
    # dump_time sweep) -- in that regime
    #   pull_func(t) - pull_func(start_time)
    #     = pull/eft * (t - start_time) * X0
    #       + pulse(t) - pulse(start_time)
    # i.e. the ramp_time/2 offset cancels.  For start_time < ramp_time a
    # piecewise restart formula would be required.
    type = ParsedFunction
    expression = 'pull_amount * ((t - start_time)/end_time_for_ramp) * x + (pulse_amplitude_bend * x + pulse_amplitude_compress) * (exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width)) - exp(-(start_time-pulse_center)*(start_time-pulse_center)/(pulse_width*pulse_width)))'
      # expression = 'pull_amount * ((t - start_time)/end_time_for_ramp) * x'

    symbol_names = 'pull_amount start_time end_time_for_ramp pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    symbol_values = '${pull_amount} ${start_time} ${end_time_for_ramp} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'
  []
[]

[BCs]
  [xfix]
    type = DirichletBC
    variable = disp_x
    boundary = '${xfix_bnd} bottom right'
    # boundary = '${xfix_bnd}'
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
    # PresetDisplacementSpatial (raccoon's spatially-aware variant of
    # MOOSE's PresetDisplacement; see src/bcs/PresetDisplacementSpatial.C).
    # On the restart side the mesh nodes sit on the dump-deformed
    # configuration, so the node's own x-coord is NOT the undeformed X0.
    # We route the function's spatial argument through coupled_x = X0
    # (the recovered AuxVariable carrying the original undeformed x of
    # each node), so the BC's spatial gradient matches the reference's
    # exactly.
    #
    # Using ADMatchedValueBC here instead would re-introduce the
    # 1/(beta*dt) accel-blowup problem: even when the recovered vel_y
    # initially matches the reference's analytical pull/eft*X0 to
    # within numerical noise, NewmarkAccelAux's TIMESTEP_END back-
    # computation amplifies that tiny mismatch by 1/(beta*dt) into the
    # accel, which then feeds into the next step's vel via NewmarkVelAux
    # and compounds.  PresetDisplacementSpatial runs the Newmark update
    # FORWARD, so the BC stays consistent regardless of tiny vel/accel
    # errors in the recovered state.
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
    # Strain-adjusted density: rho/J on the recovered configuration so the
    # inertia kernel (which integrates with use_displaced_mesh=false over the
    # restart's reference mesh -- which IS the dump-deformed mesh) integrates
    # the same total mass as the reference's continuing run.  Without this,
    # the inertia force is off by a factor of J_avg(t_dump).  Mirrors restart.i
    # in the recover_remesh_tet10_fbar test.
    type = ADStrainAdjustedDensityCustom
    strain_free_density = density
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
    order = THIRD
  []

  start_time = ${start_time}
  end_time = ${end_time}
  automatic_scaling = true
[]

[Outputs]
  print_linear_residuals = true
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
