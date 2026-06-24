//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "CoupledTimeDerivativeAux.h"

registerMooseObject("raccoonApp", CoupledTimeDerivativeAux);

InputParameters
CoupledTimeDerivativeAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Copies the time derivative of a coupled variable into an AuxVariable. "
      "Works for nodal or elemental AuxVariables (unlike TimeDerivativeAux, "
      "which is elemental-only).  Use order_derivative = SECOND to copy the "
      "second time derivative instead.");
  params.addRequiredCoupledVar(
      "v", "Variable whose (first or second) time derivative is to be copied.");
  MooseEnum order_derivative("FIRST SECOND", "FIRST");
  params.addParam<MooseEnum>(
      "order_derivative",
      order_derivative,
      "FIRST -> copy coupledDot(v); SECOND -> copy coupledDotDot(v).");
  return params;
}

CoupledTimeDerivativeAux::CoupledTimeDerivativeAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _coupled_dot(getParam<MooseEnum>("order_derivative") == "FIRST" ? coupledDot("v")
                                                                    : coupledDotDot("v"))
{
}

Real
CoupledTimeDerivativeAux::computeValue()
{
  return _coupled_dot[_qp];
}
