# Single-TET10 incompressible BENDING test (DYNAMIC variant).
#
# Same physics, BCs, mesh, X0-recovery infrastructure as pull_test.i but with
# Newmark/HHT inertia added so the system is solved as a transient elastic
# dynamic problem instead of quasi-static:
#   - hht_alpha / newmark_beta / newmark_gamma constants
#   - accel_*, vel_* AuxVariables (recovered from dump on the restart side)
#   - NewmarkAccelAux / NewmarkVelAux AuxKernels keeping the Newmark state
#     consistent with the latest disp
#   - ADInertialForce kernels (rho * a)
#   - ADDynamicStressDivergenceTensors instead of ADStressDivergenceTensors
#   - density material property added to bulk_properties
#
# Pair with pull_test_dynamic_restart.i.  The X0 fix (top-face BC driven by
# the frozen undeformed nodal x via ADMatchedValueBC) is preserved verbatim
# because Newmark auxkernels recompute vel/accel at TIMESTEP_END from
# whatever disp values the BC produced.
#
# Original docstring kept below:
# Single-TET10 incompressible BENDING test.
#
# One HEX8 element occupies the unit cube [0,1]^3.  The bottom face is pinned
# in y; front+back are clamped in z (plane strain); the top face is given a
# LINEAR-IN-X y-displacement so the element tilts/bends instead of stretching
# uniformly.  The material is Neo-Hookean (CNH) with a very small gap from
# nu=0.5 (default 0.4999) so the bulk modulus dwarfs the shear modulus --
# the standard volumetric-locking regime.
#
# Why bending, not uniaxial pull?
#   Uniform face BCs on a single TET10 give uniform F across all 8 QPs, so the
#   per-QP determinant J is constant inside the element and F-bar averaging is
#   a no-op (fbar_correction = 0 identically).  The bending BC here puts a
#   linear-in-x kinematic gradient on the top face -- top-right corner moves
#   up by `pull_amount`, top-left stays put.  The grad(u) field then varies
#   linearly across the element so F (and J) varies per QP, and F-bar actively
#   averages it.  Combined with plane-strain z and the near-incompressible
#   material, this loading drives the F-bar correction terms to large values.
#
# What this test demonstrates:
#   - With volumetric_locking_correction = true (F-bar ON), the body bends
#     freely, J stays close to 1 by construction, and fbar_correction shows
#     measurable per-QP J spread because the kinematics are heterogeneous.
#   - With volumetric_locking_correction = false (CLI override), the same
#     bending kinematics drive a much higher psie because the formulation
#     overstiffens (volumetric locking) when J variation has to be honored
#     point-by-point.  Comparing the two runs is the standard locking demo.
#
# CLI knobs:
#   raccoon-opt -i pull_test.i GlobalParams/volumetric_locking_correction=false
#   raccoon-opt -i pull_test.i pull_amount=0.5
#   raccoon-opt -i pull_test.i nu=0.4999999       # closer to incompressible
#   raccoon-opt -i pull_test.i "zfix_bnd=back"    # no plane strain (free in z)
#
# Postprocessors track J = det(F) per QP via integral and max, plus the
# F-bar diagnostics that match the recover_remesh_tet10_fbar tests so the
# numbers are directly comparable.

E = 1.0
nu = 0.4999
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 1.0

# HHT-alpha / Newmark-beta time-integration parameters.
hht_alpha = -0.25
newmark_beta = '${fparse (1-hht_alpha)^2/4}'
newmark_gamma = '${fparse 1/2-hht_alpha}'

# Loading: linear ramp of top-face y-displacement up to `pull_amount`.
end_time = 1.0
dt = 0.0125
pull_amount = 0.2

# BC ramp denominator: target_uy(top) = pull_amount * (t / end_time_for_ramp) * X0.
# Defaults to end_time so standalone runs are unchanged.  dump_time_sweep_study.sh
# overrides this so the loading rate stays constant across a sweep.
end_time_for_ramp = ${end_time}

# Gaussian BC pulse (see pull_test_dynamic.i in single_hex8_incompressible for
# the rationale).  Defaults to a 5%-of-pull_amount bump centered at t=0.05.
pulse_amplitude_bend = 0.01
pulse_amplitude_compress = 0.01
pulse_center = 0.05
pulse_width = 0.025
# Quadratic-then-linear startup window so f''(t) is continuous at t=0,
# required by PresetDisplacementSpatial (else Newmark sees a step in the
# prescribed acceleration at t=0).
ramp_time = 0.1

# Output / recovery hooks.  Plain `raccoon-opt -i pull_test.i` runs with the
# defaults below and behaves like the original standalone test.  The
# nu_sweep_study.sh script overrides `tag`, `out_dir`, and `dump_time` so the
# reference can feed a restart leg via pull_test_restart.i.
out_dir = outputs
tag = ''
# Physical time at which the displaced exodusqp recovery dump is written.
# Defaults to end_time so plain runs match the standalone test; the sweep
# script extends end_time past the recovery point and pins dump_time to the
# desired recover time.
dump_time = ${end_time}

# BC face lists.  Defaults give plane-strain bending:
#   bottom y-pinned (face),
#   left   x-pinned (face) to anchor lateral motion,
#   front + back  z-pinned  -> plane strain (this is what makes the
#                              incompressibility-driven F-bar work big),
#   top    given a LINEAR-IN-X y-displacement so the element tilts/bends.
xfix_bnd = 'left'
yfix_bnd = 'bottom'
zfix_bnd = 'front back'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  # F-bar correction.  Toggle to see locking.
  volumetric_locking_correction = true
  # Required by RecoverVariables / SolutionUserObjectQP / ComputeDeformationGradient
  # for QP-index mapping in the recovery file.  TET10_4th matches FOURTH-order
  # quadrature (11 QPs / TET10 element).
  element = TET10_4th
[]

[Problem]
  type = FEProblem
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
    # Phase-field stub used by Hencky's degradation function param even when
    # no damage is solved.  Leave at zero so degradation factor = 1.
  []
  [X0]
    # Undeformed nodal x-coordinate stored as an AuxVariable, frozen at t=0
    # via FunctionIC.  The reference's mesh IS undeformed, so this is just
    # the x of each node.  Dumped to the exodusqp recovery file so the
    # restart can recover the ORIGINAL (not deformed) x and feed it into
    # its BC -- without this the restart's spatial-coordinate BC would use
    # the post-recovery deformed mesh's coordx, which has already laterally
    # contracted from the reference's t=dump_time state.
    family = LAGRANGE
    order = SECOND
    [InitialCondition]
      type = FunctionIC
      function = 'x'
    []
  []
  [target_uy]
    # Per-node BC target for disp_y on the top face, computed from X0 so it
    # uses the ORIGINAL undeformed x rather than the deformed mesh's x.
    # MatchedValueBC below copies this into disp_y on the top boundary.
    family = LAGRANGE
    order = SECOND
  []

  # ----- Newmark / HHT dynamic state ----------------------------------------
  # accel_*, vel_* are written by NewmarkAccelAux / NewmarkVelAux at
  # TIMESTEP_END based on the latest disp_*.  The restart's epsol UO recovers
  # them via SolutionIC so the body picks up its momentum at start_time.
  [accel_x]
    order = SECOND
    family = LAGRANGE
  []
  [vel_x]
    order = SECOND
    family = LAGRANGE
  []
  [accel_y]
    order = SECOND
    family = LAGRANGE
  []
  [vel_y]
    order = SECOND
    family = LAGRANGE
  []
  [accel_z]
    order = SECOND
    family = LAGRANGE
  []
  [vel_z]
    order = SECOND
    family = LAGRANGE
  []
[]

[AuxKernels]
  [compute_target_uy]
    # Diagnostic AuxVariable (BC now consumes pull_func directly via
    # PresetDisplacementSpatial).  Kept so postprocessors that reference
    # target_uy still work.
    type = ParsedAux
    variable = target_uy
    coupled_variables = 'X0'
    constant_names       = 'pull_amount end_time_for_ramp ramp_time pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    constant_expressions = '${pull_amount} ${end_time_for_ramp} ${ramp_time} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'
    expression = 'if(t <= ramp_time, pull_amount / (2*end_time_for_ramp*ramp_time) * t*t * X0, pull_amount / end_time_for_ramp * (t - ramp_time/2) * X0) + (pulse_amplitude_bend * X0 + pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))'
    use_xyzt = true
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR'
  []

  # ----- Newmark / HHT dynamic state updates --------------------------------
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
  [x]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
    alpha = ${hht_alpha}
  []
  [y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
    alpha = ${hht_alpha}
  []
  [z]
    type = ADDynamicStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
    alpha = ${hht_alpha}
  []

  # ----- Inertia (rho * a) --------------------------------------------------
  [inertia_x]
    type = ADInertialForce
    variable = disp_x
    density = density
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
    density = density
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
    density = density
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_z
    acceleration = accel_z
  []
[]

[Functions]
  [pull_func]
    # Piecewise quadratic-then-linear ramp + Gaussian pulse.  Consumed by the
    # PresetDisplacementSpatial BC below.  See hex8/pull_test_dynamic.i for
    # the rationale (continuous f''(t) required by Newmark forward update).
    type = ParsedFunction
    expression = 'if(t <= ramp_time, pull_amount / (2*end_time_for_ramp*ramp_time) * t*t * x, pull_amount / end_time_for_ramp * (t - ramp_time/2) * x) + (pulse_amplitude_bend * x + pulse_amplitude_compress) * exp(-(t-pulse_center)*(t-pulse_center)/(pulse_width*pulse_width))'
    symbol_names = 'pull_amount end_time_for_ramp ramp_time pulse_amplitude_bend pulse_amplitude_compress pulse_center pulse_width'
    symbol_values = '${pull_amount} ${end_time_for_ramp} ${ramp_time} ${pulse_amplitude_bend} ${pulse_amplitude_compress} ${pulse_center} ${pulse_width}'
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
    preset = true
  []
  [ypull]
    # PresetDisplacementSpatial: raccoon's spatially-aware variant of MOOSE's
    # PresetDisplacement (see hex8/pull_test_dynamic.i [ypull] for the full
    # rationale -- TL;DR: ADMatchedValueBC + NewmarkAccelAux back-computing
    # a_n+1 blows up acceleration at small dt; PresetDisplacementSpatial
    # runs Newmark FORWARD instead).
    type = PresetDisplacementSpatial
    variable = disp_y
    boundary = top
    function = pull_func
    beta = ${newmark_beta}
    velocity = vel_y
    acceleration = accel_y
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
    prop_names = 'K G density'
    prop_values = '${K} ${G} ${rho}'
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

  # ----- F-bar diagnostic (same definitions as recover_remesh_tet10_fbar) ---
  # J_F     = det(Fnobar)   -- raw per-QP determinant
  # J_Fbar  = det(deformation_gradient)  -- F-bar averaged; constant per elem
  # fbar_correction          = |J/J_avg - 1|       (dimensionless kinematic spread)
  # fbar_pressure_correction = K * |J_avg - J|     (stress-scale, units of stress)
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

  end_time = ${end_time}
  automatic_scaling = true
[]

[RecoverVariables]
  [rec]
    # Both stretch_tensor (raw U) and stretch_tensor_fbar (U_bar) are dumped
    # so the restart side can exercise either approach:
    #   Approach A (NEW): consume U_bar directly
    #   Approach B (OLD): consume raw U and apply F-bar averaging at INITIAL
    # nu_sweep_study.sh overrides this list per method so only one of the two
    # stretch tensors is written per recovery file (libmesh's exodus reader
    # shadows the shorter name when both are present).
    tensor_materials = 'stress rotation_tensor stretch_tensor stretch_tensor_fbar'
    # Jacobian: per-QP J_raw at dump time, consumed by the restart's
    # adj_density_rec for exact per-QP inertia mass matching.  See
    # single_hex8_incompressible/dynamic_recovery_fix.pdf.
    materials = 'Jacobian'
  []
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = ${out_dir}/pull_test_out${tag}
    # Undisplaced output keeps the file to a single per-timestep series; with
    # use_displaced=true Exodus rewrites the mesh every step and emits a per-
    # step file (.e-s002, .e-s003, ...).  Paraview can still warp by disp_*
    # when visualizing.
    use_displaced = false
  []
  [exodusqp]
    # Displaced recovery dump consumed by pull_test_restart.i.  Written as
    # ONE slice at sync_times = dump_time so the file holds the kinematic
    # state at exactly that physical time (gated via sync_only = true).
    # use_displaced = true is required for SolutionUserObjectQP's centroid
    # lookup on the restart side; the restart's mesh comes from
    # FileMeshGenerator(file = this dump) and has to be in the deformed
    # configuration the restart will see.
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
