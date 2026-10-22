# Recover leg of the single-TET10 IMPLICIT-dynamics recovery test.
# Overlays base.i.  Mesh comes from mechorig.i's displaced dump; the kinematic
# state (R,U_bar), stress, per-QP Jacobian, and Newmark accel/vel are recovered
# through a SolutionUserObjectQP and seeded in.  Uses the HHT recover-aware force
# kernel (ADDynamicStressDivergenceTensorsRecover) so the first step's sigma_old
# comes from the recovered reference stress instead of zero.

recover_file = /home/det12/projects/raccoon/test/tests/single_tet10_implicit/tet10_imp_dump.e
out_file = /home/det12/projects/raccoon/test/tests/single_tet10_implicit/tet10_imp_rec
csv_path = /home/det12/projects/raccoon/test/tests/single_tet10_implicit/tet10_imp_rec

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
    # Approach A: consume U_bar (stretch_tensor_fbar) directly.
    tensor_materials = 'stress stretch_tensor_fbar rotation_tensor'
    materials = 'Jacobian'
    system_variables = 'd_corr X0 accel_x accel_y accel_z vel_x vel_y vel_z'
    nodal_variable_order = SECOND
    use_displaced_mesh = true
    execute_on = 'INITIAL'
    timestep = LATEST
  []
[]

[AuxVariables]
  # Reseed the recovered fields onto the (SECOND-order) targets.  order = SECOND
  # is CRITICAL so SolutionIC keeps the TET10 mid-edge node values.
  [d_corr]
    [InitialCondition]
      type = SolutionIC
      from_variable = d_corr
      variable = d_corr
      solution_uo = epsol
    []
  []
  [X0]
    [InitialCondition]
      type = SolutionIC
      from_variable = X0
      variable = X0
      solution_uo = epsol
    []
  []
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

[Kernels]
  # Swap the plain HHT stress divergence for the recover-aware variant: at
  # t_step=1 it reads sigma_old from the SolutionUserObject (recovered reference
  # stress at start_time) instead of the zero MOOSE _stress_old, avoiding the
  # missing -alpha*sigma_old force on the first restart step.
  [x]
    type = ADDynamicStressDivergenceTensorsRecover
    solution = epsol
  []
  [y]
    type = ADDynamicStressDivergenceTensorsRecover
    solution = epsol
  []
  [z]
    type = ADDynamicStressDivergenceTensorsRecover
    solution = epsol
  []
  # Inertia uses the per-QP recovery-adjusted density rho/J_recovered so the
  # mass matches the reference's per-QP mass on the deformed dump mesh.
  [inertia_x]
    density = adj_density_rec
  []
  [inertia_y]
    density = adj_density_rec
  []
  [inertia_z]
    density = adj_density_rec
  []
[]

[Functions]
  # Incremental top-y displacement = pull_func(t) - pull_func(start_time): the
  # same absolute-time profile the reference applies, minus its value at the
  # handoff, so the BC is continuous across the recovery boundary.
  [pull_func]
    type = ParsedFunction
    expression = 'pulse_amplitude*exp(-((t-pulse_center)*(t-pulse_center))/(pulse_width*pulse_width))
                - pulse_amplitude*exp(-((start_time-pulse_center)*(start_time-pulse_center))/(pulse_width*pulse_width))'
    symbol_names = 'start_time pulse_amplitude pulse_center pulse_width'
    symbol_values = '${start_time} ${pulse_amplitude} ${pulse_center} ${pulse_width}'
  []
[]

[Materials]
  # Recover F from the dumped (R, U_bar) via polar decomposition.
  [defgrad]
    recover = true
    solution = epsol
    recover_mode = polar_decomposition
  []
  # HHT recover-side stress history: stress_sol = recovered reference stress at
  # start_time; stress_old_store_sol is an identity placeholder (only used
  # multiplied by zeta = 0, so harmless).
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
  # Per-QP J_raw recovery -> adj_density_rec = rho / J_recovered for exact
  # inertia mass matching (see single_hex8_incompressible/dynamic_recovery_fix).
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
[]

[Executioner]
  start_time = ${start_time}
[]

[Outputs]
  [exodusqp]
    file_base = '${out_file}'
  []
  [exodus]
    file_base = '${out_file}_trim'
  []
[]
