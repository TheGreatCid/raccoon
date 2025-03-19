//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPenaltyConstraint.h"
#include "MooseUtils.h"
#include "RaccoonUtils.h"
#include "metaphysicl/raw_type.h"
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
  : ADKernel(parameters),
    BaseNameInterface(parameters),
    _penalty(getParam<Real>("penalty_param")),
    _u_old(_var.dofValuesOld())
{
}

ADReal
ADPenaltyConstraint::computeQpResidual()
{
  ADReal delta_d = _u[_qp] - _u_old[_qp];
  return _penalty * _test[_i][_qp] * std::pow(RaccoonUtils::Macaulay(-delta_d), 2);
}
