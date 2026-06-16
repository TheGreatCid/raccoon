# Restart leg of the F-bar comparison recover/remesh test
# (PLASTIC + QUASI-STATIC variant -- no inertia, no Newmark, no HHT).
#
# Strip relative to restart.i:
#   - Inertia kernels and ADStrainAdjustedDensityCustom removed
#   - Newmark accel / vel AuxKernels removed
#   - accel_*, vel_* AuxVariables (with SolutionIC) removed
#   - epsol drops accel_* / vel_* from system_variables
#   - Stress divergence: ADDynamicStressDivergenceTensorsRecover -> ADStressDivergenceTensors
#   - PresetDisplacement ypull -> FunctionDirichletBC tracking ypull_func_restart
#   - hht_alpha / newmark_beta / newmark_gamma constants removed
# Plasticity (J2 + JC + T-coupling) and the F-bar diagnostic survive intact.
# Pair with reference_plastic_qs.i.
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

# Plastic QS: HHT-alpha and Newmark coefficients removed (no dynamic terms).

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
    # Plastic QS: accel_* / vel_* dropped from system_variables (no Newmark).
    system_variables = 'T d d_old d_corr'
    materials = 'effective_plastic_strain'
    # Both stretch_tensor (raw U) and stretch_tensor_fbar (U_bar) are
    # available so Approach A and Approach B can both be exercised.
    tensor_materials = 'stress be_bar stretch_tensor_fbar rotation_tensor'
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
  [T]
    order = SECOND
    [InitialCondition]
      type = SolutionIC
      from_variable = T
      variable = T
      solution_uo = epsol
    []
  []
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
  # Plastic QS: accel_* and vel_* AuxVariables removed (no Newmark).
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

  # Plastic QS: Newmark accel/vel AuxKernels removed.
[]

[Kernels]
  # Plastic QS: inertia kernels removed.  Stress divergence uses the plain
  # ADStressDivergenceTensors (no HHT-alpha, no dynamic-recover first-step
  # branch).  The recovered F_bar (from defgrad/recover=true) still enters via
  # _stress through the constitutive material chain.
  [x]
    type = ADStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
    absolute_value_vector_tags = 'ref'
  []
  [y]
    type = ADStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
    absolute_value_vector_tags = 'ref'
  []
  [z]
    type = ADStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
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
    # Absolute reference top-y position vs time, same as reference_plastic_qs.i.
    type = ParsedFunction
    expression = 'if(t<=trans,v/(2*trans)*t*t,v*t-v*trans/2)'
    symbol_names = 'trans v'
    symbol_values = '${trans_time} ${final_velocity}'
  []
  [ypull_func_restart]
    # INCREMENTAL top-y displacement on the restart's deformed-config mesh.
    # The mesh comes from FileMeshGenerator(file = recover_file) -- its nodes
    # are already at the reference's positions at start_time, so disp_restart=0
    # leaves the body at that configuration.  To keep tracking the reference's
    # continued ramp the restart's BC must impose ypull_func(t) - ypull_func(s)
    # where s = start_time.  Without this the QS restart would freeze at
    # start_time while the reference keeps loading.
    type = ParsedFunction
    expression = '(if(t<=trans,v/(2*trans)*t*t,v*t-v*trans/2)) - (if(s<=trans,v/(2*trans)*s*s,v*s-v*trans/2))'
    symbol_names = 'trans v s'
    symbol_values = '${trans_time} ${final_velocity} ${start_time}'
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
  [coalescence]
    type = ADParsedMaterial
    property_name = coal
    material_property_names = 'effective_plastic_strain'
    constant_names = 'c1 c2 c3'
    constant_expressions = '${c1} ${c2} ${c3}'
    expression = ((1-c1)/(1+exp((effective_plastic_strain-c2)/c3))+c1)
    outputs = exodus
    output_properties = 'coal'
  []
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
  # Plastic QS: ADStrainAdjustedDensityCustom removed (nothing reads adj_density
  # once the inertia kernels are gone).
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
  # [JC]
  #   type = PowerLawHardening
  #   exponent = 2
  #   phase_field = d_corr
  #   reference_plastic_strain = 0.5
  #   degradation_function = nodeg
  #   yield_stress = 800
  # []
  [JC]
    type = JohnsonCookHardening
    T = T
    taylor_quinney_factor = ${Q}
    T0 = 280
    sigma_0 = 1
    reference_plastic_strain = 1
    reference_plastic_strain_rate = 1e-6
    phase_field = d_corr
    degradation_function = g
    A = 791.2
    B = 509.51
    C = 0.014
    n = 0.26
    m = 1.03
    Tm = 1033
    output_properties = 'plastic_heat_generation'
    outputs = exodus
  []
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
    plasticity_model = J2
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
    # Plastic QS: FunctionDirichletBC continuing the reference's ramp via
    # ypull_func_restart (see Functions block above for the offset rationale).
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ypull_func_restart
    preset = false
  []
[]

[Postprocessors]
  [psie_corr_active_int]
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
  # F-bar correction diagnostics (same definitions as in reference.i).
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
