# Shared base for the single-TET10 IMPLICIT-dynamics recovery test
# (Newmark/HHT).  Implicit analog of single_tet10_explicit/base.i.
#
# Differences vs the explicit (central-difference) sibling:
#   - Newmark/HHT inertia: ADDynamicStressDivergenceTensors (with alpha) plus
#     ADInertialForce (rho*a), NOT MassMatrix + ExplicitMixedOrder.  No 'mass'
#     matrix tag, no lumped mass.
#   - accel_*/vel_* are the Newmark state, updated by NewmarkAccelAux /
#     NewmarkVelAux at TIMESTEP_END; they are FULL-step quantities co-located
#     with disp, so the energy balance is time-consistent without the leapfrog
#     product-form KE the explicit test needs.
#   - Implicit NEWTON solve with nl tolerances; dt is not CFL-limited so it can
#     be much larger than the explicit dt.
#
# Overlaid by mechorig.i (reference) or mechrec.i (recover).

E = 6.9e10
nu = 0.33
K = '${fparse E/3/(1-2*nu)}'
G = '${fparse E/2/(1+nu)}'
rho = 2700

# HHT-alpha / Newmark-beta time integration.  alpha = 0 gives the classic
# average-acceleration (trapezoidal) rule: unconditionally stable and
# energy-conserving for linear response, so total_energy stays flat.
hht_alpha = 0
newmark_beta = '${fparse (1-hht_alpha)^2/4}'
newmark_gamma = '${fparse 1/2-hht_alpha}'

# Implicit -> no CFL limit, dt can be large.
dt = 4.7e-6
end_time = 1.2e-3

# Recovery handoff time.  mechorig.i dumps a synced frame here; mechrec.i
# starts from it.  Shared so both legs stay consistent.
start_time = 0.0008

# Pulse params -- used by both legs' pull_func.  The ramp params
# (pull_amount, end_time_for_ramp) live in mechorig.i since only the
# reference ramps; the recover BC is the pulse increment only.
pulse_amplitude = 0.001
pulse_center = 18e-5
pulse_width = 6e-5

xfix_bnd = 'left'
yfix_bnd = 'bottom'
zfix_bnd = 'front back'

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  volumetric_locking_correction = true
  element = TET10_4th
[]

[Problem]
  type = FEProblem
[]

# Mesh supplied by the overlay (mechorig.i = generated; mechrec.i = the dump).
[Mesh]
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
  [d]
    order = SECOND
  []
  [qual_frob]
    order = CONSTANT
    family = MONOMIAL
  []
  [frob_time]
    order = CONSTANT
    family = MONOMIAL
  []

  # Phase-field stub, kept zero; the recovery infrastructure expects it.
  [d_corr]
    order = SECOND
  []
  # Original undeformed nodal x-coordinate, frozen at t=0.  Dumped so the
  # recover can recover it; not used by the uniform BC here.
  [X0]
    family = LAGRANGE
    order = FIRST
  []

  # ----- Newmark / HHT dynamic state ---------------------------------------
  # accel_*/vel_* written by NewmarkAccelAux / NewmarkVelAux at TIMESTEP_END
  # from the latest disp_*.  The recover's epsol UO reseeds them via SolutionIC.
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

  # ----- Energy-balance diagnostics (Newmark form) -------------------------
  # kinetic_energy : 1/2 rho |v^n|^2 with v the full-step Newmark velocity.
  # psie_old       : strain energy density (kept for parity with the explicit
  #                  file; for Newmark disp and vel are co-located so
  #                  total_energy uses the current-step psie directly).
  [kinetic_energy]
    order = CONSTANT
    family = MONOMIAL
  []
  [psie_old]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # ----- Newmark accel/vel updates -----------------------------------------
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

  # ----- Energy densities ---------------------------------------------------
  [Ke]
    type = ADKineticEnergyAux
    newmark_velocity_x = vel_x
    newmark_velocity_y = vel_y
    newmark_velocity_z = vel_z
    density = density
    variable = kinetic_energy
    # INITIAL so the recover's handoff step reports the recovered KE, not 0.
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [psie_old]
    type = ADMaterialRealAux
    property = psie_active
    variable = psie_old
    execute_on = 'INITIAL TIMESTEP_END'
  []

  # ----- Mesh-quality diagnostics (parity with explicit base.i) ------------
  [qual_frob]
    type = MaterialRealAux
    property = Frobenius_norm
    variable = qual_frob
    execute_on = 'TIMESTEP_END'
  []
  [frob_time]
    type = ParsedAux
    variable = frob_time
    expression = 'if(t<0.4,0,qual_frob)'
    use_xyzt = true
    coupled_variables = qual_frob
  []
[]

[Kernels]
  # HHT-alpha stress divergence (internal force).
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

  # ----- Inertia (rho * a), Newmark -----------------------------------------
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
  # pull_func stub -- the reference (mechorig.i) and recover (mechrec.i) each
  # supply their own expression (absolute vs incremental).
  [pull_func]
    type = ParsedFunction
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
    # Newmark-forward preset displacement (raccoon's spatially-aware variant of
    # MOOSE's PresetDisplacement).  The uniform pull_func has no x-dependence,
    # so this drives the whole top face to pull_func(t) while keeping the
    # Newmark accel/vel consistent (avoids the small-dt acceleration blow-up of
    # back-computing a_n+1 from a MatchedValueBC).
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
    # density is AD here (ADInertialForce / ADKineticEnergyAux read an AD rho).
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
[]

[Postprocessors]
  [frob]
    type = ElementExtremeValue
    variable = frob_time
    value_type = max
    execute_on = 'TIMESTEP_END'
  []

  # ---- Energy balance: total_energy = psi(u^n) + 1/2 rho |v^n|^2 -----------
  # For Newmark, disp and vel are co-located at the full step, so both terms
  # are taken at the current step (no leapfrog dt/2 offset).  Displaced-mesh
  # integrals with constant density keep them continuous across the recovery.
  [psie_active_int]
    type = ADElementIntegralMaterialProperty
    mat_prop = psie_active
    use_displaced_mesh = true
    #   execute_on = 'INITIAL TIMESTEP_END'
  []
  [kinetic_energy_int]
    type = ElementIntegralVariablePostprocessor
    variable = kinetic_energy
    use_displaced_mesh = true
    #   execute_on = 'INITIAL TIMESTEP_END'
  []
  [total_energy]
    type = ParsedPostprocessor
    expression = 'psie_active_int + kinetic_energy_int'
    pp_names = 'psie_active_int kinetic_energy_int'
    #   execute_on = 'INITIAL TIMESTEP_END'
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
  [disp_y_top_int]
    type = SideAverageValue
    variable = disp_y
    boundary = top
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
    tensor_materials = 'stress rotation_tensor stretch_tensor stretch_tensor_fbar'
    materials = 'Jacobian'
  []
[]

[Outputs]
  print_linear_residuals = false
  time_step_interval = 1
  [exodusqp]
    type = Exodus
    use_displaced = false
  []
  [exodus]
    type = Exodus
  []
  [csv]
    type = CSV
    file_base = ${csv_path}
  []
[]
