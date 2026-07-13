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


# CRITICAL: must match the reference's end_time_for_ramp.  The incremental
# BC formula `pull_amount * (t - start_time) / end_time_for_ramp` is
# derived from `ref_total(t) - ref_total(start_time)`, which only cancels
# correctly if both sides use the same end_time_for_ramp.  Reference uses
# end_time_for_ramp = ${end_time} = 1; restart must too.


[Mesh]
    file = ${recover_file}
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
    system_variables = 'd X0 accel_x accel_y accel_z vel_x vel_y vel_z'
    use_displaced_mesh = true
    execute_on = 'INITIAL'
    timestep = LATEST
  []
[]


[AuxVariables]
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
  [d_corr]
    [InitialCondition]
      type = SolutionIC
      from_variable = d
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

[]

[AuxKernels]

  [ke_lf]
    density = rho_const
  []
[]

[Kernels]

  [mass_x]
    density = adj_density
  []
  [mass_y]
    density = adj_density
  []
  [mass_z]
    density = adj_density
  []
[]

[Functions]
 
  [ypull_func_restart]
    # INCREMENTAL top-y displacement = pull_func(t) - pull_func(start_time),
    # i.e. the SAME absolute-time profile the reference applies (mechorig.i's
    # pull_func) minus its value at the handoff.  The recover mesh sits at the
    # reference's start_time position, so subtracting pull_func(start_time)
    # keeps the BC continuous across the recovery boundary while reproducing
    # the reference's subsequent motion exactly.
    type = ParsedFunction
    expression = '(pulse_amplitude_compress)*exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))
                + pull_amount*min(1,max(0,(t-pulse_center-4*pulse_width)/end_time_for_ramp))
                - (pulse_amplitude_compress)*exp(-(start_time-pulse_center)*(start_time-pulse_center)/(pulse_width*pulse_width))
                - pull_amount*min(1,max(0,(start_time-pulse_center-4*pulse_width)/end_time_for_ramp))'
    symbol_names = 'ramp_time pull_amount start_time end_time_for_ramp pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    symbol_values = '${ramp_time} ${pull_amount} ${start_time} ${end_time_for_ramp} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'
  []
[]

[BCs]
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
    recover_apply_fbar_to_U = false
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
[]

[UserObjects]
  [Terminator]
    expression = '1<0'
  []
[]

[Executioner]

 # [TimeIntegrator]
 #   type = ExplicitMixedOrder
 #   mass_matrix_tag = 'mass'
 #   second_order_vars = 'disp_x disp_y disp_z'
 #   use_constant_mass = true
 # []
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
#    type = ConstantDT
#    dt = ${dt}
  []

 start_time = ${start_time}
  end_time = ${end_time}
[]
[Outputs]
  [exodusqp]
        type = Exodus
        file_base = '${out_file}'
  []
  [exodus]
        type = Exodus
        file_base = '${out_file}_trim'
  []
#       [exfail]
#         type = Exodus
#         file_base = '/gpfs/detorre/exodusfiles/zerotriax/tet10s/epd/ztq_thick_fail'
#         execute_on = 'FAILED'
#   []  
[]

