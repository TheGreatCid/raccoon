//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "LargeDeformationJ2PlasticModified.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", LargeDeformationJ2PlasticityModified);

InputParameters
LargeDeformationJ2PlasticityModified::validParams()
{
  InputParameters params = LargeDeformationPlasticityModel::validParams();
  params.addClassDescription("Large deformation $J_2$ plasticity. The exponential constitutive "
                             "update is used to update the plastic deformation.");
  params.addRequiredCoupledVar("T", "The temperature of the domain");
  params.addParam<MaterialPropertyName>(
      "softening_param_del", "delta", "Thermal softneing parameter");
  params.addParam<MaterialPropertyName>("thermal_conduct", "kappa", "Thermal conductivity");
  params.addParam<MaterialPropertyName>("ref_temp", "T0", "Reference temperature");
  params.addParam<MaterialPropertyName>("ref_pstrain", "ep0", "Reference plastic strain");
  params.addParam<MaterialPropertyName>(
      "ref_pstrain_rate", "ep0_dot", "Reference plastic strain rate");
  params.addParam<MaterialPropertyName>("ref_stress", "sigma0", "Reference Stress");
  params.addParam<MaterialPropertyName>("softening_param_k", "k", "Softening parameter");
  return params;
}

LargeDeformationJ2PlasticityModified::LargeDeformationJ2PlasticityModified(
    const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    _T(adCoupledValue("T")),
    _del(getParam<Real>("softening_param_del")),
    _kappa(getParam<Real>("thermal_conduct")),
    _T0(getParam<Real>("ref_temp")),
    _ep0(getParam<Real>("ref_pstrain")),
    _ep0_dot(getParam<Real>("ref_pstrain_rate")),
    _sigma0(getParam<Real>("ref_stress")),
    _k(getParam<Real>("softening_param_k"))
{
}

void
LargeDeformationJ2PlasticityModified::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
{
  // First assume no plastic increment
  ADReal delta_ep = 0;
  Fe = Fe * _Fp_old[_qp].inverse();
  stress = _elasticity_model->computeMandelStress(Fe);

  // Compute the flow direction following the Prandtl-Reuss flow rule.
  // We guard against zero denominator.
  ADRankTwoTensor stress_dev = stress.deviatoric();
  ADReal stress_dev_norm = stress_dev.doubleContraction(stress_dev);
  if (MooseUtils::absoluteFuzzyEqual(stress_dev_norm, 0))
    stress_dev_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
  stress_dev_norm = std::sqrt(1.5 * stress_dev_norm);
  _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;

  // Return mapping
  ADReal phi = computeResidual(stress_dev_norm, delta_ep);
  if (phi > 0)
    returnMappingSolve(stress_dev_norm, delta_ep, _console);
  _ep[_qp] = _ep_old[_qp] + delta_ep;
  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  // Update stress
  Fe = Fe * delta_Fp.inverse();
  stress = _elasticity_model->computeCauchyStress(Fe);
  _hardening_model->plasticEnergy(_ep[_qp]);
}

Real
LargeDeformationJ2PlasticityModified::computeReferenceResidual(
    const ADReal & effective_trial_stress, const ADReal & delta_ep)
{
  return raw_value(
      effective_trial_stress -
      _elasticity_model->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
          .doubleContraction(_Np[_qp]));
}

ADReal
ADgdeg_T::rexp_calc()
{
  return 1 - _del * (std::exp((_T[_qp] - _T0) / _kappa) - 1);
}

ADReal
ADgdeg_T::rtanh_calc()
{
  Real kab = std::log(2) / (_kappa * (std::log((1 + _del) / _del)) * (1 + _del) - 1);
  return 1 - std::tanh(kab * (_T[_qp] - _T0));
}

ADReal
LargeDeformationJ2PlasticityModified::computeResidual(const ADReal & effective_trial_stress,
                                                      const ADReal & delta_ep)
{
  if (rexp_calc() > rtanh_calc())
  {
    ADReal r = ADgdeg_T::rexp_calc();
  }
  else
  {
    ADReal r = ADgdeg_T::rtanh_calc();
  }
  return (1 / _dt) * (_ep[_qp] - ep_old[_qp]) -
         std::pow(elasticity_model
                          ->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
                          .doubleContraction(_Np[_qp]) /
                      (_sigma0 * std::pow(1 + _ep[_qp] / _ep0, 2) * _r),
                  2);
}

ADReal
LargeDeformationJ2PlasticityModified::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                        const ADReal & delta_ep)
{
  return -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
              .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 2);
}
