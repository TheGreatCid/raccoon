[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = ../gold/kal.m
  []
[]

[Variables]
  [temp]
  []
[]

[Kernels]
  [conduction]
    type = HeatConduction
    variable = 'temp'
  []
  [del]
    type = ADCoefMatSource
    variable = 'temp'
    prop_names = 'del_delt'
  []
[]

[temp_surround]
  type = DirichletBC
  variable = 'temp'
  boundary = 'bottom load other'
  value = '293' #K
[]

[Materials]
  [del_delt]
    type = ADParsedMaterial

  []

[]
