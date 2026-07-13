# Restart leg of the single-HEX8 explicit-dynamics (central-difference)
# recovery test.  Pairs with pull_test_dynamic.i (reference).
#
# The restart's mesh is the reference's exodus dump (FileMeshGenerator);
# disp_x/y/z start at 0 (incremental).  ExplicitMixedOrderRecover seeds the
# integrator's solutionUDot / solutionUDotDot from the dump's nodal vel_* /
# accel_* AuxVariables, so the central-difference solve picks up at the
# reference's state at dump_time instead of starting at rest.
#
# F at INITIAL is seeded from (R, U_bar) via ComputeDeformationGradient
# (recover = true), same pattern as the implicit Newmark variant.

E = 1.0
nu = 0.3
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 1.0

dt = 0.001
start_time = 0.5
end_time = 1
pull_amount = 1   # must match the reference's pull_amount
# CRITICAL: must match the reference's end_time_for_ramp.  The incremental
# BC formula `pull_amount * (t - start_time) / end_time_for_ramp` is
# derived from `ref_total(t) - ref_total(start_time)`, which only cancels
# correctly if both sides use the same end_time_for_ramp.  Reference uses
# end_time_for_ramp = ${end_time} = 1; restart must too.
end_time_for_ramp = 1

out_dir = outputs
tag = ''
output_tag = ${tag}
recover_file = ${out_dir}/pull_test_out_disp${tag}.e
# recover_timestep = 0.5

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
    # Approach C reads raw stretch_tensor (U from polar decomp of F_dump_raw)
    # and recomputes F-bar on the new mesh from F_total = F_inc * F_dump_raw.
    tensor_materials = 'stretch_tensor_fbar stress stretch_tensor rotation_tensor'
    # Per-QP Jacobian recovered so adj_density_rec (= rho / J_recovered) gets
    # the correct per-QP mass.  See dynamic_recovery_fix.pdf for why this
    # matters: det(_F_NoFbar) in approach A is element-constant J_bar, not
    # the reference's actual per-QP J_raw.
    materials = 'Jacobian'
    # CD recovery: ExplicitMixedOrderRecover reads vel_*/accel_* from this
    # SolutionUserObject via directValue(node, var_name) and writes into the
    # integrator's UDot/UDotDot vectors.  X0 / d_corr recovered as well so
    # the same dump file works for future spatial-BC restart inputs.
    system_variables = 'd_corr X0 accel_x accel_y accel_z vel_x vel_y vel_z'
    use_displaced_mesh = true
    execute_on = 'INITIAL'
    timestep = LATEST
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
    [InitialCondition]
      type = SolutionIC
      from_variable = d_corr
      variable = d_corr
      solution_uo = epsol
    []
  []
  [X0]
    family = LAGRANGE
    order = FIRST
    [InitialCondition]
      type = SolutionIC
      from_variable = X0
      variable = X0
      solution_uo = epsol
    []
  []

  # vel_*/accel_* on the restart side carry the CD integrator state mirrored
  # back out for downstream postprocessing / chained recovery dumps.  They
  # don't drive the integrator -- ExplicitMixedOrderRecover seeds the
  # integrator's UDot/UDotDot directly from the dump via directValue.  These
  # AuxVars get repopulated each step via TimeDerivativeAux.
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

  # Product-form (leapfrog-conserved) kinetic energy density,
  # 1/2 rho (v^{n+1/2} . v^{n-1/2}); see LeapfrogKineticEnergyAux.
  [kinetic_energy_lf]
    order = CONSTANT
    family = MONOMIAL
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

  # Leapfrog product-form kinetic energy: 1/2 rho0 (v^{n+1/2} . v^{n-1/2}).  The
  # aux reads the OLD state of vel_* for v^{n-1/2}.  density = rho_const (a plain
  # CONSTANT rho0, NOT adj_density): using the same constant density as the base
  # run makes the displaced-mesh KE integral a physical integral over Omega(t),
  # continuous across the handoff.  (adj_density = rho0/J(0.5) would scale the
  # restart KE by a constant J(0.5) relative to the base -- the mismatch we saw.)
  [Ke_lf]
    type = LeapfrogKineticEnergyAux
    velocity_x = vel_x
    velocity_y = vel_y
    velocity_z = vel_z
    density = rho_const
    variable = kinetic_energy_lf
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
    # INCREMENTAL top-y displacement past start_time.  Reference's BC was
    # pull_amount * (t/end_time_for_ramp); restart's mesh sits at the
    # reference's t=start_time position, so the BC must add the increment.
    type = ParsedFunction
    expression = 'pull_amount * ((t - start_time) / end_time_for_ramp)'
    symbol_names = 'pull_amount start_time end_time_for_ramp'
    symbol_values = '${pull_amount} ${start_time} ${end_time_for_ramp}'
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
    # Same caveat as the reference: no spatial BC variant yet; the function
    # is uniform across the top face.
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
    # Approach C: recover raw U from the dump, compute F_total = F_inc * F_dump_raw
    # at every step, apply F-bar to F_total ONCE with the corrected J_avg formula
    # (using recovered J_dump_raw and change-of-variables to the original volume).
    # Required for explicit dynamics -- approach A/B's multiplicative F-bar
    # composition gives ~1e-4 rel_err that Newton can absorb in implicit but CD
    # cannot.  Mutually exclusive with recover_apply_fbar_to_U.
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
    # Non-AD constant rho, consumed by adj_density_rec below.  MassMatrix
    # itself reads adj_density_rec.
    type = GenericConstantMaterial
    prop_names = 'rho_const'
    prop_values = '${rho}'
  []
  # Use the GEOMETRIC det(_F_store_noFbar) (= the dump-deformed mesh's local
  # Jacobian via polar decomp of the recovered raw F) for density adjustment.
  # That's the same J MOOSE uses for _JxW, so (rho/J) * _JxW cancels at
  # floating-point level, not just analytically.  Bypasses the J_recovered
  # (from constitutive _F.det()) vs J_mesh_geometric mismatch that was
  # causing the 5e-4 lumped-mass discrepancy.
  # StrainAdjustedDensityCustom with base_name='adj' declares 'adj_density'.
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
    # Should equal ref's psie_active_int at t=dump_time (= start_time here)
    # if the kinematic recovery is bit-exact.
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
    execute_on = 'INITIAL'
  []

  # ---- Leapfrog energy balance (product-form KE + strain energy at u^n) ----
  # H^n = 1/2 rho0 (v^{n+1/2}.v^{n-1/2}) + psi(u^n): the product-form KE is
  # centered at step n, so it pairs with the strain energy at the start-of-step
  # displacement u^n, evaluated at TIMESTEP_BEGIN (material and displaced mesh
  # both at u^n).  Plain DISPLACED-mesh integrals with a CONSTANT density
  # (rho_const on the KE) -> physical integrals over Omega(t), continuous with
  # the base run across the handoff.  NOTE: the first restart step has
  # vel_*_old = 0 (no history across the recover), so KE_lf is a one-step
  # transient there -- ignore it.
  [kinetic_energy_lf_int]
    type = ElementIntegralVariablePostprocessor
    variable = kinetic_energy_lf
    use_displaced_mesh = true
  []
  [psie_old_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [total_energy_lf]
    type = ParsedPostprocessor
    expression = 'psie_old_int + kinetic_energy_lf_int'
    pp_names = 'psie_old_int kinetic_energy_lf_int'
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
  # SideAverageValue removed -- it triggered creation of a boundary CDG
  # (defgrad_face) which then errored on missing stretch_tensor_xx_n
  # (boundary materials have vlc forced to false, taking the need_U_raw
  # branch which the dump doesn't have).  If a side-integrated diagnostic
  # is needed later, add it after the recovery is verified end-to-end.
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
    # Skip mass matrix recomputation each step -- it's constant for this test
    # (density, mesh, and adj_density_rec don't change).  Eliminates a
    # potential source of step-to-step numerical noise.
    use_constant_mass = true
  []

  [TimeStepper]
    type = ConstantDT
    dt = ${dt}
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
