# Writes a synthetic phase-field crack band on a TET10 box, as ../planar_crack/band.i does, with a
# crack-opening displacement that stretches the band across the crack line y = 0.5 (x < x_tip):
# d = exp(-dist / l), and disp_y jumps smoothly by 2 * u0 across the line.

a = 0
x_tip = 0.6
l = 0.1
amp = 1
u0 = 0.05

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 12
    ny = 12
    nz = 1
    zmax = 0.1
    elem_type = TET10
  []
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [d]
    order = SECOND
  []
  [disp_x]
    order = SECOND
  []
  [disp_y]
    order = SECOND
  []
  [disp_z]
    order = SECOND
  []
  [cx]
    order = SECOND
  []
  [cy]
    order = SECOND
  []
  [cz]
    order = SECOND
  []
  [e]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [d]
    type = FunctionIC
    variable = d
    function = band
  []
  [disp_x]
    type = FunctionIC
    variable = disp_x
    function = '${fparse 0.5 * u0} * y'
  []
  [disp_y]
    type = FunctionIC
    variable = disp_y
    function = opening
  []
  [cx]
    type = FunctionIC
    variable = cx
    function = x
  []
  [cy]
    type = FunctionIC
    variable = cy
    function = y
  []
  [cz]
    type = FunctionIC
    variable = cz
    function = z
  []
  [e]
    type = FunctionIC
    variable = e
    function = 'x + 2 * y + 3 * z'
  []
[]

[Functions]
  [band]
    type = ParsedFunction
    expression = 'yc := 0.5 + a * sin(pi * min(x, x_tip));
                  dist := if(x < x_tip, abs(y - yc), sqrt((x - x_tip)^2 + (y - yc)^2));
                  min(1, amp * exp(-dist / l))'
    symbol_names = 'a x_tip l amp'
    symbol_values = '${a} ${x_tip} ${l} ${amp}'
  []
  [opening]
    type = ParsedFunction
    expression = 'yc := 0.5 + a * sin(pi * min(x, x_tip));
                  u0 * tanh((y - yc) / l) * if(x < x_tip, 1, exp(-(x - x_tip) / l))'
    symbol_names = 'a x_tip l u0'
    symbol_values = '${a} ${x_tip} ${l} ${u0}'
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  [exodus]
    type = Exodus
    file_base = band
  []
[]
