//* This file is part of the RACCOON framework

#include "RotateBC.h"

registerMooseObject("raccoonApp", RotateBC);

InputParameters
RotateBC::validParams()
{
  InputParameters params = ADDirichletBCBase::validParams();
  params.addRequiredCoupledVar("theta", "Variable providing angles at nodal values");
  params.addClassDescription("Fixes rotation in the xy plane");
  params.addRequiredParam<MooseEnum>("component", MooseEnum("x y", "x"), "x or y component");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");

  return params;
}

RotateBC::RotateBC(const InputParameters & parameters)
  : ADDirichletBCBase(parameters),
    _theta(coupledValueOld("theta")),
    _xy(getParam<MooseEnum>("component"))
{
}

ADReal
RotateBC::computeQpValue()
{
  const VariableValue * x_disp = &coupledDofValues("displacements", 0);
  const VariableValue * y_disp = &coupledDofValues("displacements", 1);
  if (_xy == "x")
  {
    return (*x_disp)[_qp] / tan(_theta[_qp]);
  }
  if (_xy == "y")
  {
    return (*y_disp)[_qp] * tan(_theta[_qp]);
  }
  return 0;
}
