//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADTimeDerivativeCustom.h"

registerMooseObject("raccoonApp", ADTimeDerivativeCustom);

InputParameters
ADTimeDerivativeCustom::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription("outputting the time derivative");
  params.addRequiredCoupledVar("deriv_var", "variable to take time derivative of");

  return params;
}

ADTimeDerivativeCustom::ADTimeDerivativeCustom(const InputParameters & parameters)
  : ADKernel(parameters), _d_dot(adCoupledDot("deriv_var"))
{
}

ADReal
ADTimeDerivativeCustom::computeQpResidual()
{
  return _test[_i][_qp] * (_u[_qp] - _d_dot[_qp]);
}
