# Reference leg of the single-TET10 IMPLICIT-dynamics recovery test.
# Overlays base.i.  Newmark/HHT transient elastic dynamics; dumps a synced
# displaced frame at start_time for mechrec.i to recover from.

out_file = /home/det12/projects/raccoon/test/tests/single_tet10_implicit/tet10_imp_dump
csv_path = /home/det12/projects/raccoon/test/tests/single_tet10_implicit/tet10_imp_orig

# Ramp loading params -- reference-only (the recover BC is pulse-increment only).
pull_amount = 0.0025
end_time_for_ramp = 0.6

[Mesh]
  file = 'mesh_box.e'
  type = FileMesh
  construct_node_list_from_side_list = False
  # [gmg]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 2
  #   ny = 2
  #   nz = 2
  #   xmin = 0
  #   xmax = 1
  #   ymin = 0
  #   ymax = 1
  #   zmin = 0
  #   zmax = 1
  #   elem_type = TET10
  # []
[]

[AuxVariables]
  [X0]
    family = LAGRANGE
    order = FIRST
    [InitialCondition]
      type = FunctionIC
      function = 'x'
    []
  []
[]

[Functions]
  [pull_func]
    # Absolute (uniform) top-face y-displacement: linear ramp to pull_amount
    # plus a Gaussian pulse.  No spatial dependence -- the whole top face moves
    # together (uniaxial, matching the explicit test's ExplicitFunctionDirichletBC).
    type = ParsedFunction
    expression = 'pull_amount*(t/end_time_for_ramp)
                + pulse_amplitude*exp(-((t-pulse_center)*(t-pulse_center))/(pulse_width*pulse_width))'
    symbol_names = 'pull_amount end_time_for_ramp pulse_amplitude pulse_center pulse_width'
    symbol_values = '${pull_amount} ${end_time_for_ramp} ${pulse_amplitude} ${pulse_center} ${pulse_width}'
  []
[]

[Outputs]
  # Recovery dump: a single DISPLACED-mesh frame at start_time, so mechrec.i can
  # FileMesh-load it as the recover run's deformed reference.
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
