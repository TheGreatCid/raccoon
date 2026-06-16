# Reference simulation for the F-bar comparison recover/remesh test
# (ELASTIC-ONLY variant).
#
# Strip relative to reference.i:
#   - No plasticity (no [J2], no [JC] hardening), no plastic_heat_generation
#   - No temperature (T variable, heat kernels, JC's T coupling) -- all gone
#   - No `coalescence` material (depended on effective_plastic_strain)
#   - RecoverVariables drops `be_bar` and `effective_plastic_strain`
#   - Postprocessors drop psip_active_int / ep_int
# Pair with restart_elastic.i.  Outputs land in a separate directory so the
# elastic and plastic sweeps don't clobber each other.
#
# Original reference comment kept below:
# Reference simulation for the F-bar comparison recover/remesh test.
#
# Tuned to maximize F-bar correction work so the recover-side choice between
# (A) recovering U_bar directly and (B) recovering raw U + re-correcting on the
# new mesh produces a measurable difference in step-1 state.
#
# Tunings vs recover_remesh_tet10/reference.i:
#   - Poisson ratio nu = 0.49           (near-incompressible -> high K/G ratio
#                                        -> volumetric locking severe -> F-bar
#                                        correction does heavy lifting)
#   - sigma_0 = 0.1 (was 1)             (lower yield -> earlier and stronger
#                                        plastic flow -> isochoric plastic
#                                        deformation drives det F variation
#                                        between QPs that F-bar must average)
#   - final_velocity = 0.5 (was 0.05)   (10x deformation magnitude)
#   - end_time = 2.0 (was 1)            (more steps to accumulate F-bar history)
#   - The two output exports include `Fnobar` and `deformation_gradient`
#     (= F-bar) so the per-QP F vs F-bar gap is visible in Paraview.
#
# Otherwise mirrors recover_remesh_tet10/reference.i: TET10, FOURTH-order
# quadrature (11 QPs/elem), J2 plasticity with Johnson-Cook hardening,
# temperature, and a phase-field fracture sub-app.

E = 201.8e3
nu = 0.49
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 7900
Gc = 6
l = 0.2
psic = 15
Q = 0.9
specific_heat = 4.47e-4
thermal_conductivity = 4.4e-4
c1 = 0.1
c2 = 0.9
c3 = 0.1

Tinit = 293

hht_alpha = -0.25
newmark_beta = '${fparse (1-hht_alpha)^2/4}'
newmark_gamma = '${fparse 1/2-hht_alpha}'

end_time = 1.5
dt = 0.1

# Mesh-refinement knob (overridable from CLI: raccoon-opt -i reference.i n=4).
# Controls all three element-counts: nx = ny = nz = ${n}.
n = 1

# All Exodus / CSV outputs land in this subdirectory.  Overridable from CLI:
#   raccoon-opt -i reference.i out_dir=runs/mesh4
out_dir = outputs

# Suffix appended to file_base for every Output block.  Sweeps that vary a
# parameter at fixed n (e.g. nu, final_velocity) use this to label each run's
# files so they don't overwrite each other.  Empty by default.  CLI:
#   raccoon-opt -i reference.i tag=_nu_0p49
tag = ''

# Physical time at which the displaced exodusqp recovery dump is written.
# Defaults to end_time so the FINAL dump matches the original behaviour.
# Convergence scripts that extend end_time past ref_end_time override this
# to point at the original ref_end_time so the recovery dump stays at the
# physically-correct deformation state.  CLI:
#   raccoon-opt -i reference.i end_time=2.0 dump_time=1.5
dump_time = ${end_time}

# BC constraint pattern (CLI-overridable for bc_sweep_study.sh).
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
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${n}
    ny = ${n}
    nz = ${n}
    xmin = 0
    xmax = 3
    ymin = 0
    ymax = 3
    zmin = 0
    zmax = 3
    elem_type = TET10
  []
[]

# [MultiApps]
#   [fracture]
#     type = TransientMultiApp
#     input_files = 'fracture_ref.i'
#     cli_args = 'Gc=${Gc};l=${l};out_file=reference_fbar_d'
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
  #   initial_condition = ${Tinit}
  # []
[]

[AuxVariables]
  # Elastic variant: no T AuxVariable -- nothing reads it.
  [d]
    order = SECOND
  []
  [d_old]
    order = SECOND
  []
  [d_corr]
    order = SECOND
  []

  [accel_x]
    order = SECOND
  []
  [vel_x]
    order = SECOND
  []
  [accel_y]
    order = SECOND
  []
  [vel_y]
    order = SECOND
  []
  [accel_z]
    order = SECOND
  []
  [vel_z]
    order = SECOND
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
    density = density
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
    density = density
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
    density = density
    use_displaced_mesh = false
    beta = ${newmark_beta}
    gamma = ${newmark_gamma}
    alpha = ${hht_alpha}
    velocity = vel_z
    acceleration = accel_z
    absolute_value_vector_tags = 'ref'
  []
  [x]
    type = ADDynamicStressDivergenceTensors
    variable = disp_x
    component = 0
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    absolute_value_vector_tags = 'ref'
  []
  [y]
    type = ADDynamicStressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = true
    alpha = ${hht_alpha}
    absolute_value_vector_tags = 'ref'
  []
  [z]
    type = ADDynamicStressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = true
    alpha = ${hht_alpha}
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

[RecoverVariables]
  [rec]
    # Elastic variant: be_bar dropped (no plastic flow), effective_plastic_strain
    # dropped (no plasticity).  stretch_tensor / stretch_tensor_fbar both written
    # so the restart side can exercise Approach A and Approach B.
    tensor_materials = 'stress rotation_tensor stretch_tensor stretch_tensor_fbar'
    materials = ''
  []
[]

[Materials]
  # Elastic variant: coalescence removed (depended on effective_plastic_strain).
  [defgrad]
    type = ComputeDeformationGradient

    output_properties = 'deformation_gradient Fnobar'
    outputs = exodus
  []
  [bulk_properties]
    type = ADGenericConstantMaterial
    prop_names = 'K G l Gc psic density thermal_conductivity specific_heat'
    prop_values = '${K} ${G} ${l} ${Gc} ${psic} ${rho} ${thermal_conductivity} ${specific_heat}'
  []
  [reg_density]
    type = MaterialADConverter
    ad_props_in = 'density'
    reg_props_out = 'reg_density'
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
  # Elastic variant: [J2] and [JC] (plasticity + hardening) removed.
  # [stress] consumes only the elasticity model.
  [stress]
    type = ComputeLargeDeformationStress
    elasticity_model = hencky
  []

  # ----- F-bar diagnostic ---------------------------------------------------
  # F-bar replaces F at each QP with (J_avg/J)^(1/3) * F, where J = det(F) and
  # J_avg is the per-element mean.  Thus:
  #   det(Fnobar)               = J         (per-QP raw determinant)
  #   det(deformation_gradient) = J_avg     (constant within an element)
  # The pointwise correction magnitude is |J/J_avg - 1| = |J_F/J_Fbar - 1|.
  # If F-bar is doing no work it is identically zero; large values mean the
  # volumetric-locking correction is active.
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
  # Stress-scale F-bar correction.  At each QP this is the volumetric Cauchy
  # stress DIFFERENCE between (a) the F-bar response, which uses J_avg, and
  # (b) what the same constitutive law would return on the unmodified J:
  #     fbar_pressure_correction ~ K * (J_avg - J)
  # Unlike fbar_correction (a dimensionless kinematic spread), this carries
  # the bulk modulus K and therefore GROWS with nu in the way you'd expect
  # of an "amount of volumetric correction" diagnostic.
  [fbar_pressure_correction]
    type = ADParsedMaterial
    property_name = fbar_pressure_correction
    material_property_names = 'J_F J_Fbar K'
    expression = K*abs(J_Fbar-J_F)
    outputs = exodus
  []
[]

trans_time = 1.0
final_velocity = 0.2

[Functions]
  [ypull_func]
    type = ParsedFunction
    expression = 'if(t<=trans,v/(2*trans)*t*t,v*t-v*trans/2)'
    symbol_names = 'trans v'
    symbol_values = '${trans_time} ${final_velocity}'
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
    function = ypull_func
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
  # F-bar correction diagnostics.
  # fbar_correction_int    - integrated |J/J_avg - 1| over Omega_0; scales with
  #                          how much volumetric averaging F-bar is doing.
  # fbar_correction_max    - worst-case per-QP correction magnitude anywhere
  #                          in the domain.
  # J_F_int, J_Fbar_int    - sanity check: integrals of det(F) and det(F-bar)
  #                          over Omega_0 must agree by construction of F-bar.
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

[Dampers]
  [jac]
    type = ElementJacobianDamper
    max_increment = 0.1
  []
[]

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

  fixed_point_max_its = 25
  fixed_point_rel_tol = 1e-7
  fixed_point_abs_tol = 1e-8
  accept_on_max_fixed_point_iteration = true
  abort_on_solve_fail = false
[]

[Outputs]
  print_linear_residuals = false
  [exodus]
    type = Exodus
    file_base = ${out_dir}/reference_fbar_out_${n}${tag}
    use_displaced = false
  []
  [exodusqp]
    type = Exodus
    file_base = ${out_dir}/reference_fbar_out_disp_${n}${tag}
    # use_displaced = true is required for SolutionUserObjectQP's
    # centroid-to-centroid lookup on the restart side -- the restart's mesh
    # comes from FileMeshGenerator(file = this dump), so the dump's node
    # coordinates have to be in the deformed configuration the restart will
    # see.  Exodus stores ONE mesh definition per file, so this output only
    # writes ONCE -- gated to a single sync time = dump_time so the dump
    # captures the kinematic state at exactly that physical time.  Default
    # dump_time = end_time matches the original FINAL behaviour; convergence
    # scripts extend end_time past ref_end_time and override
    # dump_time = ref_end_time so the restart still recovers from the
    # original ref_end_time slice.
    use_displaced = true
    sync_times = '${dump_time}'
    sync_only = true
  []
  [csv]
    type = CSV
    file_base = ${out_dir}/reference_fbar_out_${n}${tag}
  []
[]
