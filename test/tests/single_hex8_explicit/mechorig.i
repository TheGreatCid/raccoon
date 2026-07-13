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


[Mesh]
    file = 'mesh.e'
##  [shift]
#    type = TransformGenerator
#    input = fmg
#    transform = TRANSLATE
#    vector_value = '0 0 0.2'
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
    #expression = 'if(t <= ramp_time, pull_amount / (2*end_time_for_ramp*ramp_time) * t*t * x, pull_amount / end_time_for_ramp * (t - ramp_time/2) * x) + (pulse_amplitude_bend * x + pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))'
   
   
   # expression = '(pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))'
    expression = '(pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))+pull_amount * min(1,max(0,(t-pulse_center-4*pulse_width)/end_time_for_ramp))'
    symbol_names = 'pull_amount end_time_for_ramp ramp_time pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    symbol_values = '${pull_amount} ${end_time_for_ramp} ${ramp_time} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'
  []
 # [pull_func]
 #   # Linear ramp.  No spatial dependence (top face moves uniformly).
 #   type = ParsedFunction
 #   expression = 'pull_amount * (t / end_time_for_ramp)'
 #   symbol_names = 'pull_amount end_time_for_ramp'
 #   symbol_values = '${pull_amount} ${end_time_for_ramp}'
 # []
[]

[BCs]
  [ypull]
    boundary = top

    function = pull_func
  []
[]



[Executioner]

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

[Outputs]
  # Recovery dump: a single DISPLACED-mesh slice at start_time (= dump_time),
  # so base.i + mechrec.i can FileMesh-load it as the recover run's deformed
  # reference (mirrors the proven pull_test_dynamic.i dump pattern).
  [exodusqp]
        file_base = '${out_file}'
        use_displaced = true
        sync_times = '${start_time}'
        sync_only = true
  []
  [exodus]
    file_base = '${out_file}_trim'
  []

[]

