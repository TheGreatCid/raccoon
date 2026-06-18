# Restart leg of the F-bar comparison recover/remesh test
# (ELASTIC-ONLY variant).
#
# Strip relative to restart.i:
#   - No plasticity, no hardening, no JC/T coupling
#   - No `coalescence` material
#   - epsol drops `be_bar` from tensor_materials and drops `materials` line
#   - [stress] uses only elasticity_model
#   - psip_active_int and ep_int postprocessors removed
# Pair with reference_elastic.i.
#
# Original restart comment kept below:
# Restart leg of the F-bar comparison recover/remesh test.
#
# Reads reference.i's output via SolutionUserObjectQP on the same TET10 mesh.
# Reads BOTH stretch_tensor (= raw U) and stretch_tensor_fbar (= U-bar) so
# whichever path ComputeDeformationGradient is configured to consume,
# the recovery data is available.
#
# To compare Approach A (recover U-bar directly) vs Approach B (recover raw U,
# re-derive F-bar on the new mesh): rerun this input with the corresponding
# branch / parameter setting in ComputeDeformationGradient and compare the
# postprocessor CSV at start_time + dt.  The reference run was tuned so that
# F-bar correction does meaningful work (nu = 0.49, sigma_0 = 0.1,
# final_velocity = 0.5).

E = 201.8e3
nu = 0.49
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 7900
Gc = 6
l = 0.2
psic = 15
Q = 0
specific_heat = 4.47e-4
thermal_conductivity = 4.4e-4
c1 = 0.1
c2 = 0.9
c3 = 0.1

hht_alpha = -0.25
newmark_beta = '${fparse (1-hht_alpha)^2/4}'
newmark_gamma = '${fparse 1/2-hht_alpha}'

# Must match reference.i
# ref_end_time = 1.5
dt = 0.1
start_time = 1.5
end_time = 2.5

trans_time = 1.0
final_velocity = 0.5

# Mesh-refinement knob: must match the n used in reference.i so the recover
# filename and SolutionUserObjectQP mesh line up.
# Overridable from CLI:  raccoon-opt -i restart.i n=4
n = 1

# All Exodus / CSV outputs land in this subdirectory.  Must match the out_dir
# used by reference.i so the recover_file path resolves.
out_dir = outputs

# Suffix used to locate the reference's recovery dump.  Must match the `tag`
# passed to the reference.i run that produced the dump.  CLI:
#   raccoon-opt -i restart.i tag=_nu_0p49
tag = ''

# Suffix appended to file_base for every Output block written by this restart.
# Defaults to `tag` so a single-method run keeps shared naming with the
# reference; sweeps that exercise multiple recovery methods (e.g. approach A
# vs approach B) against the same reference dump override this to a more
# specific value like '_nu_0p49_old' / '_nu_0p49_new' so each method's
# outputs don't clobber each other.  CLI:
#   raccoon-opt -i restart.i tag=_nu_0p49 output_tag=_nu_0p49_old
output_tag = ${tag}

recover_file = ${out_dir}/reference_fbar_out_disp_${n}${tag}.e

# Which time slice of recover_file to read into SolutionUserObjectQP.  LATEST
# is fine when reference and restart line up at the end of the reference run,
# but the convergence script needs to point at the slice corresponding to
# ref_end_time precisely (so it can extend the reference one dt past
# ref_end_time and still recover from the correct state).  CLI override:
#   raccoon-opt -i restart.i recover_timestep=16
recover_timestep = LATEST

# BC constraint pattern (CLI-overridable for bc_sweep_study.sh).  Must match
# the values passed to the reference run.
xfix_bnd = 'left top'
yfix_bnd = 'bottom'
zfix_bnd = 'front back'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
  element = TET10_4th
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
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
    # Elastic variant: T dropped from system_variables.
    system_variables = 'd d_old d_corr accel_x accel_y accel_z vel_x vel_y vel_z'
    # Elastic variant: `materials = 'effective_plastic_strain'` removed
    # (no plasticity) and `be_bar` dropped from `tensor_materials`.
    # Jacobian: per-QP J_raw for adj_density_rec (see dynamic_recovery_fix.pdf).
    materials = 'Jacobian'
    tensor_materials = 'stress stretch_tensor_fbar rotation_tensor'
    nodal_variable_order = SECOND
    use_displaced_mesh = true
    execute_on = 'INITIAL'
    timestep = ${recover_timestep}
  []
[]

# [MultiApps]
#   [fracture]
#     type = TransientMultiApp
#     input_files = 'fracture_restart.i'
#     cli_args = 'Gc=${Gc};l=${l};start_time=${start_time};recover_file=${recover_file};out_file=restart_fbar_d'
#     execute_on = 'TIMESTEP_END'
#     clone_parent_mesh = no
#   []
# []

# [Transfers]
#   [to_coal]
#     type = MultiAppCopyTransfer
#     variable = coal
#     source_variable = coal
#     to_multi_app = fracture
#   []
#   [from_d]
#     type = MultiAppCopyTransfer
#     variable = d
#     source_variable = d
#     from_multi_app = fracture
#   []
#   [to_psie]
#     type = MultiAppCopyTransfer
#     variable = psie_corr_active
#     source_variable = psie_corr_active
#     to_multi_app = fracture
#   []
#   [to_psip]
#     type = MultiAppCopyTransfer
#     variable = psip_active
#     source_variable = psip_active
#     to_multi_app = fracture
#   []
# []

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
  # [T]
  #   order = SECOND
  #   [InitialCondition]
  #     type = SolutionIC
  #     from_variable = T
  #     variable = T
  #     solution_uo = epsol
  #   []
  # []
[]

[AuxVariables]
  # Elastic variant: T AuxVariable dropped (nothing reads it).
  # [d]
  #   order = SECOND
  #   [InitialCondition]
  #     type = SolutionIC
  #     from_variable = d
  #     variable = d
  #     solution_uo = epsol
  #   []
  # []
  # [d_old]
  #   order = SECOND
  #   [InitialCondition]
  #     type = SolutionIC
  #     from_variable = d_old
  #     variable = d_old
  #     solution_uo = epsol
  #   []
  # []
  [d_corr]
    order = SECOND
    # [InitialCondition]
    #   type = SolutionIC
    #   from_variable = d_corr
    #   variable = d_corr
    #   solution_uo = epsol
    # []
  []
  [accel_x]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_x
      variable = accel_x
      solution_uo = epsol
    []
  []
  [vel_x]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_x
      variable = vel_x
      solution_uo = epsol
    []
  []
  [accel_y]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_y
      variable = accel_y
      solution_uo = epsol
    []
  []
  [vel_y]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_y
      variable = vel_y
      solution_uo = epsol
    []
  []
  [accel_z]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = accel_z
      variable = accel_z
      solution_uo = epsol
    []
  []
  [vel_z]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = vel_z
      variable = vel_z
      solution_uo = epsol
    []
  []
[]

[AuxKernels]
  # [d_old]
  #   type = CopyValueAux
  #   source = d
  #   variable = d_old
  #   execute_on = 'TIMESTEP_END'
  # []
  # [d_corr]
  #   type = ParsedAux
  #   variable = d_corr
  #   coupled_variables = 'd d_old'
  #   expression = 'min(1,max(d_old,max(0,d)))'
  # []

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
    absolute_value_vector_tags = 'ref'
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
    absolute_value_vector_tags = 'ref'
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
    absolute_value_vector_tags = 'ref'
  []
  [x]
    type = ADDynamicStressDivergenceTensorsRecover
    variable = disp_x
    component = 0
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    solution = epsol
    absolute_value_vector_tags = 'ref'
  []
  [y]
    type = ADDynamicStressDivergenceTensorsRecover
    variable = disp_y
    component = 1
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    solution = epsol
    absolute_value_vector_tags = 'ref'
  []
  [z]
    type = ADDynamicStressDivergenceTensorsRecover
    variable = disp_z
    component = 2
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    solution = epsol
    absolute_value_vector_tags = 'ref'
  []
  # [hcond_time]
  #   type = ADHeatConductionTimeDerivative
  #   variable = T
  #   density_name = density
  #   specific_heat = specific_heat
  #   absolute_value_vector_tags = 'ref'
  # []
  # [hcond]
  #   type = ADHeatConduction
  #   variable = T
  #   thermal_conductivity = thermal_conductivity
  #   absolute_value_vector_tags = 'ref'
  # []
  # [heat_source]
  #   type = ADCoefMatSource
  #   variable = T
  #   coefficient = -1
  #   prop_names = 'plastic_heat_generation'
  #   absolute_value_vector_tags = 'ref'
  # []
[]

[Functions]
  [ypull_func]
    type = ParsedFunction
    expression = 'if(t<=trans,v/(2*trans)*t*t,v*t-v*trans/2)'
    symbol_names = 'trans v'
    symbol_values = '${trans_time} ${final_velocity}'
  []
[]

[Materials]
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
  # Elastic variant: coalescence removed (depended on effective_plastic_strain).
  [defgrad]
    type = ComputeDeformationGradient
    recover = true
    solution = epsol
    recover_mode = polar_decomposition
    # Higham iterative polar decomposition -- tighter precision than the default
    # eigendecomposition.  Trying as a precision-sensitivity test for whether
    # the OLD recovery method (approach B) converges with sharper R.
    use_iterative_polar_decomposition = true
    output_properties = 'deformation_gradient Fnobar'
    outputs = exodus
  []
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc psic density thermal_conductivity specific_heat'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic} ${rho} ${thermal_conductivity} ${specific_heat}'
  []
  [dens]
    # Kept for diagnostics (det(_F_NoFbar) which in approach A is element-
    # constant J_bar).  Inertia consumes adj_density_rec below, which uses
    # the per-QP J_raw recovered from the reference's dump.
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
  [g]
    type = PowerDegradationFunction
    phase_field = d_corr
    property_name = g
    expression = (1-d_corr)^p*(1-eta)+eta
    parameter_names = 'p eta'
    parameter_values = '2 1e-6'
  []
  [hencky]
    type = CNHIsotropicElasticity
    bulk_modulus = K
    shear_modulus = G
    phase_field = d_corr
    degradation_function = g
    decomposition = NONE
  []
  # Elastic variant: [J2] and [JC] removed.  [stress] uses elasticity only.
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
  []

  # ----- F-bar diagnostic ---------------------------------------------------
  # Same definitions as in reference.i: pointwise |J/J_avg - 1| where
  #   J     = det(Fnobar)               (raw per-QP determinant)
  #   J_avg = det(deformation_gradient) (constant on the element by F-bar)
  # Comparing fbar_correction_int / _max between restart and reference at the
  # same physical time tells you whether the recovered + re-corrected F-bar
  # reproduces the same level of volumetric averaging on the restart side.
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
  # Stress-scale F-bar correction (carries K, so it grows with nu):
  #     fbar_pressure_correction ~ K * |J_avg - J|
  # See reference.i for the rationale.
  [fbar_pressure_correction]
    type = ADParsedMaterial
    property_name = fbar_pressure_correction
    material_property_names = 'J_F J_Fbar K'
    expression = K*abs(J_Fbar-J_F)
    outputs = exodus
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
    type = PresetDisplacement
    variable = disp_y
    boundary = top
    function = '0'
    beta = ${newmark_beta}
    velocity = vel_y
    acceleration = accel_y
  []
[]

[Postprocessors]
  [psie_corr_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
  []
  # Elastic variant: psip_active_int and ep_int postprocessors removed.
  # F-bar correction diagnostics (same definitions as in reference_elastic.i).
  [fbar_correction_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = fbar_correction
    use_displaced_mesh = true
  []
  [fbar_correction_max]
    type = ADElementExtremeMaterialProperty
    mat_prop = fbar_correction
    value_type = max
  []
  [fbar_pressure_correction_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = fbar_pressure_correction
    use_displaced_mesh = true
  []
  [fbar_pressure_correction_max]
    type = ADElementExtremeMaterialProperty
    mat_prop = fbar_pressure_correction
    value_type = max
  []
  [J_F_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = J_F
    use_displaced_mesh = true
  []
  [J_Fbar_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = J_Fbar
    use_displaced_mesh = true
  []
[]

# [Dampers]
#   [jac]
#     type = ElementJacobianDamper
#     max_increment = 0.1
#   []
# []

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'

  [TimeIntegrator]
    type = ImplicitEuler
  []

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 100

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

  abort_on_solve_fail = false
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = ${out_dir}/restart_fbar_out_${n}${output_tag}
    use_displaced = false
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [csv]
    type = CSV
    file_base = ${out_dir}/restart_fbar_out_${n}${output_tag}
    # TIMESTEP_END only.  EXEC_INITIAL would write a row at t=start_time but
    # MOOSE calls Material::initStatefulProperties (not the regular
    # computeQpProperties) at INITIAL, so material-property integrals are
    # zero at that flag for any material that doesn't override the stateful
    # init.  The mesh-convergence script therefore compares the restart's
    # FIRST TIMESTEP_END row (t = start_time + dt) to a reference row at
    # the same time -- which requires the reference to be extended by one
    # dt past ref_end_time (see mesh_convergence.sh).
    execute_on = 'TIMESTEP_END'
  []
[]
