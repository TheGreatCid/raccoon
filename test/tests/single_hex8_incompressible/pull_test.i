# Single-HEX8 incompressible BENDING test.
#
# One HEX8 element occupies the unit cube [0,1]^3.  The bottom face is pinned
# in y; front+back are clamped in z (plane strain); the top face is given a
# LINEAR-IN-X y-displacement so the element tilts/bends instead of stretching
# uniformly.  The material is Neo-Hookean (CNH) with a very small gap from
# nu=0.5 (default 0.4999) so the bulk modulus dwarfs the shear modulus --
# the standard volumetric-locking regime.
#
# Why bending, not uniaxial pull?
#   Uniform face BCs on a single HEX8 give uniform F across all 8 QPs, so the
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

# Loading: linear ramp of top-face y-displacement up to `pull_amount`.
end_time = 1.0
dt = 0.0125
pull_amount = 0.2

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

# BC ramp denominator: target_uy(top) = pull_amount * (t / end_time_for_ramp) * X0.
# Defaults to end_time so standalone runs are unchanged.  dump_time_sweep_study.sh
# overrides this so the loading rate stays constant across a sweep that varies
# the dump time.
end_time_for_ramp = ${end_time}

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
  # for QP-index mapping in the recovery file.  HEX8_3rd matches THIRD-order
  # quadrature (8 QPs / HEX8 element).
  element = HEX8_3rd
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
  [d_corr]
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
    order = FIRST
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
    order = FIRST
  []
[]

[AuxKernels]
  [compute_target_uy]
    type = ParsedAux
    variable = target_uy
    coupled_variables = 'X0'
    constant_names       = 'pull_amount end_time_for_ramp'
    constant_expressions = '${pull_amount} ${end_time_for_ramp}'
    expression = 'pull_amount * (t/end_time_for_ramp) * X0'
    use_xyzt = true
    # Run before the solve so the BC sees the correct target each step.
    execute_on = 'INITIAL TIMESTEP_BEGIN LINEAR'
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
[]

[Functions]
  [pull_func]
    # u_y(top) = pull_amount * (t/end_time) * x  -- LINEAR IN X.
    # At t = end_time:
    #   top-left  node (x=0): u_y = 0           (stays at y=1)
    #   top-right node (x=1): u_y = pull_amount (rises to y=1+pull_amount)
    # The interior grad(u) field then varies linearly with position, so F
    # and J vary per QP -- this is what gives F-bar non-trivial work to do.
    type = ParsedFunction
    expression = 'pull_amount * (t/end_time_for_ramp) * x'
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
    # MatchedValueBC enforces disp_y = target_uy on the top face.  target_uy
    # is itself an AuxVariable computed via ParsedAux from X0 (the original
    # undeformed x) -- so the BC's spatial gradient comes from the original
    # x, not from the deformed mesh's coordx.  This makes the restart's BC
    # behave the same as the reference's even after recovery, since X0
    # gets recovered through SolutionIC and matches the reference's value.
    type = ADMatchedValueBC
    variable = disp_y
    boundary = top
    v = target_uy
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
    prop_names = 'K G'
    prop_values = '${K} ${G}'
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
    order = SECOND
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
    materials = ''
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
