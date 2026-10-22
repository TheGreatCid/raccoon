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

csv_path = /home/det12/projects/raccoon/test/tests/single_tet10_explicit/orig

# Generated single-block TET10 mesh (ignore the missing mesh_box.e for now).
# Rename the generated sidesets to the *_ns names base.i's DirichletBCs use;
# 'top' is kept as-is for the pull BC.  The reference dump carries these
# boundary names into the exodus so the recover mesh (FileMesh of the dump)
# resolves the same BCs.
[Mesh]
  file = 'mesh_box.e'
  type = FileMesh
  construct_node_list_from_side_list = False
  #  [gmg]
  #    type = GeneratedMeshGenerator
  #    dim = 3
  #    nx = 2
  #    ny = 2
  #    nz = 2
  #    elem_type = TET10
  #  []
  #  [rename]
  #    type = RenameBoundaryGenerator
  #    input = gmg
  #    old_boundary = 'left bottom front back'
  #    new_boundary = 'left_ns bottom_ns front_ns back_ns'
  #  []
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
[]

[AuxKernels]
[]

[Kernels]
[]
[Functions]
  [pull_func]
    # Piecewise quadratic-then-linear ramp + Gaussian pulse, consumed by the
    # PresetDisplacementSpatial BC.  See hex8/pull_test_dynamic.i.
    # expression = 'if(t <= ramp_time, pull_amount / (2*end_time_for_ramp*ramp_time) * t*t * x, pull_amount / end_time_for_ramp * (t - ramp_time/2) * x) + (pulse_amplitude_bend * x + pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))'
    #    expression = 'if(t<0.32,pull_amount * (t/2) + ( pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width)),pull_amount*(0.32/2)+(pulse_amplitude_compress)*exp(-(0.32-pulse_center)*(0.32-pulse_center)/(pulse_width*pulse_width)))'
    #    expression = 'pull_amount * (t/end_time_for_ramp)+(pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))'
    #expression = 'pull_amount*(t/end_time_for_ramp)+pulse_amplitude * exp(-((t-pulse_center)*(t-pulse_center))/(pulse_width*pulse_width))'
    # symbol_names = 'pulse_amplitude pull_amount end_time_for_ramp ramp_time pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    expression = 'pulse_amplitude * exp(-((t-pulse_center)*(t-pulse_center))/(pulse_width*pulse_width))'
    # symbol_values = '${pulse_amplitude} ${pull_amount} ${end_time_for_ramp} ${ramp_time} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'

    symbol_names = 'pulse_amplitude pull_amount end_time_for_ramp pulse_center pulse_width'
    symbol_values = '${pulse_amplitude} ${pull_amount} ${end_time_for_ramp} ${pulse_center} ${pulse_width}'
  []
[]

[BCs]
  [ypull]
    boundary = top

    function = pull_func
  []
[]

[Executioner]

  [TimeIntegrator]
    type = ExplicitMixedOrderHRZ
    mass_matrix_tag = 'mass'
    second_order_vars = 'disp_x disp_y disp_z'
    use_constant_mass = true
    mass_lumping_type = hrz
  []

  [TimeStepper]
    type = ConstantDT
    dt = ${dt}
  []

  end_time = ${end_time}
  # No nl_*_tol / nl_max_its -- ExplicitMixedOrder bypasses the nonlinear solver.
[]

[Outputs]
  [exodusqp]
    file_base = '/home/det12/projects/raccoon/test/tests/single_tet10_explicit/orig'
    sync_times = '${start_time}'
    sync_only = true
  []
  [exodus]
    file_base = '/home/det12/projects/raccoon/test/tests/single_tet10_explicit/orig_trim'
  []
[]

