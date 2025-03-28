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
  params.addParam<Real>("epsilon", 1e-6, "The penalty param");
  params.addParam<bool>("smooth", false, "Whether or not to use a smoothed macaulay bracket");

  return params;
}

ADPenaltyConstraint::ADPenaltyConstraint(const InputParameters & parameters)
  : ADKernel(parameters),
    BaseNameInterface(parameters),
    _penalty(getParam<Real>("penalty_param")),
    _epsilon(getParam<Real>("epsilon")),
    _u_old(valueOld()),
    _smooth(getParam<bool>("smooth"))
{
}

ADReal
ADPenaltyConstraint::computeQpResidual()
{
  ADReal function = 0;
  ADReal delta_d = _u_old[_qp] > 0 ? _u_old[_qp] - _u[_qp] : 0 - _u[_qp];

  if (_smooth)
    function = 0.5 * (std::sqrt(delta_d * delta_d + _epsilon * _epsilon) + delta_d);
  else
    RaccoonUtils::Macaulay(delta_d);

  return -_penalty * _test[_i][_qp] * function;
}
