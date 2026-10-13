//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPFFDiffusionRecover.h"

registerMooseObject("raccoonApp", ADPFFDiffusionRecover);

InputParameters
ADPFFDiffusionRecover::validParams()
{
  InputParameters params = ADPFFDiffusion::validParams();
  params.addClassDescription("The diffusion term in the phase-field evolution equation. The weak "
                             "form is $(\\grad w, \\dfrac{2\\Gc l}{c_0} \\grad d)$.");

  // Not read by the kernel; still accepted so existing input files that set it keep working.
  params.addParam<Real>("qp", 8, "Unused.");
  params.addParam<MaterialPropertyName>("Fbarprev", "Fbar from previous simulation");
  return params;
}

ADPFFDiffusionRecover::ADPFFDiffusionRecover(const InputParameters & parameters)
  : ADPFFDiffusion(parameters), _Fbarprev(getADMaterialProperty<RankTwoTensor>("Fbarprev"))
{
}

ADReal
ADPFFDiffusionRecover::computeQpResidual()
{
  const auto mult = _grad_u[_qp] * _Fbarprev[_qp];
  const ADReal value = _grad_test[_i][_qp] * _Fbarprev[_qp] * mult;
  return 2 * _Gc[_qp] * _l[_qp] / _c0[_qp] * value;
}
