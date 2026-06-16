# Restart leg of the single-HEX8 incompressible bending recovery test
# (PLASTIC_QS variant -- QS + J2 + PowerLawHardening).
# Pairs with pull_test_plastic_qs.i (reference).
#
# Adds to the elastic QS restart (pull_test_restart.i):
#   - be_bar tensor recovered through epsol so the isochoric plastic state
#     carries over.
#   - LargeDeformationJ2PlasticityBeBar(recover = true, solution = epsol) and
#     effective_plastic_strain in epsol's materials list.
#   - PowerLawHardening called `JC` (parity with the recover_remesh tet10 test).
# No Newmark / inertia (QS).
#
# Original docstring (from the elastic QS restart) below:
# Restart leg of the single-HEX8 incompressible bending recovery test.
# Pairs with pull_test.i (reference).
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

# PowerLawHardening parameters (must match pull_test_plastic_qs.i).
yield_stress = 0.03
hardening_exponent = 2
hardening_ref_pstrain = 0.5

# Time stepping: must match the reference.i values.
dt = 0.0125
start_time = 0.5
end_time = 1.0
pull_amount = 0.2
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
    tensor_materials = 'stress be_bar stretch_tensor_fbar rotation_tensor'
    materials = 'effective_plastic_strain'
    # Recover X0 (frozen undeformed nodal x-coords) so the BC can use the
    # ORIGINAL x rather than the deformed-mesh coordx, and d_corr (phase
    # field stub).
    system_variables = 'd_corr X0'
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
[]

[AuxKernels]
  [compute_target_uy]
    type = ParsedAux
    variable = target_uy
    coupled_variables = 'X0'
    constant_names       = 'pull_amount start_time end_time_for_ramp'
    constant_expressions = '${pull_amount} ${start_time} ${end_time_for_ramp}'
    expression = 'pull_amount * X0 * (t - start_time) / end_time_for_ramp'
    use_xyzt = true
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
    # The mesh's top nodes already sit at the reference's positions at
    # start_time, so the BC has to apply ypull_func(t) - ypull_func(start_time).
    # Without this the QS restart would freeze at start_time while the
    # reference keeps ramping past it.
    type = ParsedFunction
    expression = 'pull_amount * x * ((t - start_time) / end_time_for_ramp)'
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
    # MatchedValueBC enforces disp_y = target_uy on the top face.  target_uy
    # is computed via ParsedAux from RECOVERED X0 -- i.e. the original
    # undeformed x of each node, which equals the reference's BC spatial
    # argument.  This is the X0 fix: without it the BC used the deformed
    # mesh's coordx and the top-right shortfell by ~1% of the increment.
    type = ADMatchedValueBC
    variable = disp_y
    boundary = top
    v = target_uy
  []
[]

[Materials]
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
  [JC]
    type = PowerLawHardening
    phase_field = d_corr
    degradation_function = nodeg
    yield_stress = ${yield_stress}
    exponent = ${hardening_exponent}
    reference_plastic_strain = ${hardening_ref_pstrain}
  []
  [J2]
    type = LargeDeformationJ2PlasticityBeBar
    phase_field = d_corr
    hardening_model = JC
    relative_tolerance = 1e-08
    recover = true
    solution = epsol
    output_properties = 'effective_plastic_strain'
    outputs = exodus
  []
  [coalescence]
    type = ADParsedMaterial
    property_name = coal
    material_property_names = 'effective_plastic_strain'
    constant_names = 'c1 c2 c3'
    constant_expressions = '0.1 0.9 0.1'
    expression = ((1-c1)/(1+exp((effective_plastic_strain-c2)/c3))+c1)
    outputs = exodus
    output_properties = 'coal'
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
    plasticity_model = J2
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
  [psip_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psip_active
    use_displaced_mesh = true
  []
  [ep_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = effective_plastic_strain
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
