//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPenaltyConstraint.h"
#include "MooseUtils.h"
registerMooseObject("raccoonApp", ADPenaltyConstraint);

InputParameters
ADPenaltyConstraint::validParams()
{
  InputParameters params = ADKernel::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("The diffusion term in the phase-field evolution equation. The weak "
                             "form is $(\\grad w, \\dfrac{2\\Gc l}{c_0} \\grad d)$.");

  params.addParam<Real>("penalty_param", 1e6, "The penalty param");
  return params;
}

ADPenaltyConstraint::ADPenaltyConstraint(const InputParameters & parameters)
  : ADKernel(parameters), BaseNameInterface(parameters), _penalty(getParam<Real>("penalty_param"))
{
}

ADReal
ADPenaltyConstraint::computeQpResidual()
{
  ADReal d = 0;
  if (-_u[_qp] >= 0)
    d = -(_u[_qp] * _u[_qp]);

  return _penalty * _test[_i][_qp] * d;
}
