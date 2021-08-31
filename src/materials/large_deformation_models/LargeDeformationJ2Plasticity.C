//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "LargeDeformationJ2Plasticity.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", LargeDeformationJ2Plasticity);

InputParameters
LargeDeformationJ2Plasticity::validParams()
{
  InputParameters params = LargeDeformationPlasticityModel::validParams();
  params.addClassDescription("Large deformation $J_2$ plasticity. The exponential constitutive "
                             "update is used to update the plastic deformation.");
  params.addRequiredParam<MaterialPropertyName>("ref_yield_stress",
                                                "The reference yield stress $\\sigma_0$");
  params.addRequiredParam<MaterialPropertyName>("exponent",
                                                "The exponent n in the power law hardening $n$");
  params.addRequiredParam<MaterialPropertyName>(
      "reference_plastic_strain", "The $\\epsilon_0$ parameter in the power law hardening");
  params.addRequiredParam<MaterialPropertyName>("Taylor_Quinney", "The Taylor_Quinney factor");
  params.addRequiredParam<MaterialPropertyName>("c_v", "c_v");
  params.addRequiredParam<MaterialPropertyName>("rho", "the density");
  params.addRequiredParam<MaterialPropertyName>("R", "TS_Rate");
  params.addRequiredParam<MaterialPropertyName>("T0", "ref_temp");

  return params;
}

LargeDeformationJ2Plasticity::LargeDeformationJ2Plasticity(const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    _T(declareProperty<Real>(prependBaseName("Temp"))),
    _T_old(getMaterialPropertyOldByName<Real>(prependBaseName("Temp"))),
    _T0(getADMaterialProperty<Real>(prependBaseName("T0"))),
    _sigma_0(getADMaterialProperty<Real>(prependBaseName("ref_yield_stress"))),
    _n(getADMaterialProperty<Real>(prependBaseName("exponent", true))),
    _ep0(getADMaterialProperty<Real>(prependBaseName("reference_plastic_strain", true))),
    _xi(getADMaterialProperty<Real>(prependBaseName("Taylor_Quinney"))),
    _cv(getADMaterialProperty<Real>(prependBaseName("c_v"))),
    _rho(getADMaterialProperty<Real>(prependBaseName("rho"))),
    _R(getADMaterialProperty<Real>(prependBaseName("R")))

{
}

void
LargeDeformationJ2Plasticity::initQpStatefulProperties()
{
  LargeDeformationPlasticityModel::initQpStatefulProperties();
  //_sigma_y[_qp] = _sigma_0[_qp];
  _T[_qp] = 293;
}

void
LargeDeformationJ2Plasticity::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
{
  // First assume no plastic increment
  ADReal delta_ep = 0;
  //------------------------

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

  stress_dev = stress.deviatoric();
  // std::cout << stress_dev << std::endl;
  // // Update temp
  // std::cout << std::sqrt(1.5 * (stress_dev.doubleContraction(stress_dev))) << std::endl;
  computeTemperature(
      delta_ep,
      std::sqrt(1.5 * (stress_dev.doubleContraction(stress_dev)))); //----------------
}

Real
LargeDeformationJ2Plasticity::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                       const ADReal & delta_ep)
{
  return raw_value(
      effective_trial_stress -
      _elasticity_model->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
          .doubleContraction(_Np[_qp]));
}

ADReal
LargeDeformationJ2Plasticity::computeResidual(const ADReal & effective_trial_stress,
                                              const ADReal & delta_ep)
{
  // ADReal f = _ep_old[_qp] + delta_ep;
  // ADReal S = computeSigyDeriv(delta_ep, effective_trial_stress, 0);
  // ADReal dS = computeSigyDeriv(delta_ep, effective_trial_stress, 1);
  //  ADReal ddS = computeSigyDeriv(delta_ep, effective_trial_stress, 2);
  return effective_trial_stress - // Adding chain rule term due to the temp dependence on sigy
         _elasticity_model->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
             .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 1);
  //  +
  // _ep0[_qp] * _n[_qp] * dS * (std::pow(f / _ep0[_qp] + 1, 1 / _n[_qp] + 1) - 1) *
  //     (1 / (_n[_qp] + 1));
}

ADReal
LargeDeformationJ2Plasticity::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  // Set everything as variables and rewrite to make this easier to debug
  // ADReal f = _ep_old[_qp] + delta_ep;

  // ADReal S = computeSigyDeriv(delta_ep, effective_trial_stress, 0);
  // ADReal dS = computeSigyDeriv(delta_ep, effective_trial_stress, 1);
  // ADReal ddS = computeSigyDeriv(delta_ep, effective_trial_stress, 2);
  return -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
              .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 2);

  // 2 * dS * std::pow((_sigma_0[_qp] + f) / _sigma_0[_qp], 1 / _n[_qp]) -
  // (_sigma_0[_qp] * _n[_qp] * ddS *
  //  (1 - std::pow((_sigma_0[_qp] + f) / _sigma_0[_qp], 1 / _n[_qp] + 1))) /
  //     (_n[_qp] + 1);
}

void
LargeDeformationJ2Plasticity::computeTemperature(const ADReal & delta_ep,
                                                 const ADReal & effective_stress)
{
  // std::cout << "----------------------------------------" << std::endl;

  _T[_qp] = raw_value((0.5 * effective_stress * delta_ep) / (_cv[_qp] * _rho[_qp]) + _T_old[_qp]);
  // std::cout << "eff_stress * del_ep = " << effective_stress * delta_ep << std::endl;
  // std::cout << "T_old = " << _T_old[_qp] << std::endl;
  // std::cout << raw_value(_T[_qp]) << std::endl;
}

ADReal
LargeDeformationJ2Plasticity::computeSigyDeriv(const ADReal & delta_ep,
                                               const ADReal & effective_trial_stress,
                                               const unsigned int derivative)
{
  if (derivative == 0)
  {
    return _sigma_0[_qp] * std::exp((293 - _T[_qp]) / 10);
  }

  if (derivative == 1)
  {
    return -((_sigma_0[_qp] * _ep0[_qp] * _xi[_qp] *
              std::exp((1 / _R[_qp]) *
                           (-1 * (effective_trial_stress * _xi[_qp] * (_ep_old[_qp] + delta_ep)) *
                            (1 / (_cv[_qp] * _rho[_qp]))) -
                       _T_old[_qp] + _T0[_qp])) /
             (_cv[_qp] * _rho[_qp] * _R[_qp]));
  }

  if (derivative == 2)
  {
    return ((_sigma_0[_qp] * std::pow(_ep0[_qp], 2) * std::pow(_xi[_qp], 2) *
             std::exp((1 / _R[_qp]) *
                          (-1 * (effective_trial_stress * _xi[_qp] * (_ep_old[_qp] + delta_ep)) *
                           (1 / (_cv[_qp] * _rho[_qp]))) -
                      _T_old[_qp] + _T0[_qp])) /
            (std::pow(_cv[_qp], 2) * std::pow(_rho[_qp], 2) * std::pow(_R[_qp], 2)));
  }
  mooseError(name(), "internal error: unsupported derivative order.");
  return 0;
}
