//* David's development code

#include "ADgdeg_T.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", ADgdeg_T);

InputParameters
ADgdeg_T::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription("Computes stiffness degration for a temperature solver");
  params.addParam<MaterialPropertyName>("Taylor_Quinney", "Q", "The Taylor_Quinney Factor");
  params.addParam<MaterialPropertyName>(
      "softening_param_del", "delta", "Thermal softneing parameter");
  params.addParam<MaterialPropertyName>("thermal_conduct", "kappa", "Thermal conductivity");
  params.addParam<MaterialPropertyName>("ref_temp", "T0", "Reference temperature");
  params.addParam<MaterialPropertyName>("ref_pstrain", "ep0", "Reference plastic strain");
  params.addParam<MaterialPropertyName>(
      "ref_pstrain_rate", "ep0_dot", "Reference plastic strain rate");
  params.addParam<MaterialPropertyName>("ref_stress", "sigma0", "Reference Stress");
  params.addParam<MaterialPropertyName>("effective_stress", "sigma_eff", "Effective_stress");
  params.addParam<MaterialPropertyName>("softening_param_k", "k", "Softening parameter");
  return params;
}

ADgdeg_T::ADgdeg_T(const InputParameters & parameters)
  : ADKernelValue(parameters),
    _Q(getParam<Real>("Taylor_Quinney")),
    _tau_bar(getADMaterialPropertyByName<Real>("effective_kirchoff_stress")),
    _ep(getADMaterialPropertyByName<Real>("effective_plastic_strain")),
    _del(getParam<Real>("softening_param_del")),
    _kappa(getParam<Real>("thermal_conduct")),
    _T0(getParam<Real>("ref_temp")),
    _ep0(getParam<Real>("ref_pstrain")),
    _ep0_dot(getParam<Real>("ref_pstrain_rate")),
    _sigma0(getParam<Real>("ref_stress")),
    _sigma_eff(getParam<Real>("effective_stress")),
    _k(getParam<Real>("softening_param_k"))
{
}

ADReal
ADgdeg_T::rexp_calc()
{
  return 1 - _del * (std::exp((_u[_qp] - _T0) / _kappa) - 1);
}

ADReal
ADgdeg_T::rtanh_calc()
{
  Real kab = std::log(2) / (_kappa * (std::log((1 + _del) / _del)) * (1 + _del) - 1);
  return 1 - std::tanh(kab * (_u[_qp] - _T0));
}

ADReal
ADgdeg_T::precomputeQpResidual()
{

  if (rexp_calc() > rtanh_calc())
  {
    ADReal r = ADgdeg_T::rexp_calc();
  }
  else
  {
    ADReal r = ADgdeg_T::rtanh_calc();
  }
  return _Q * _tau_bar[_qp] * _ep0 *
         std::pow(_sigma_eff / (_sigma0 * std::pow(1 + _ep[_qp] / _ep0, 2) * _r), 2);
}
