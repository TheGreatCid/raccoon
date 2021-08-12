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
  return params;
}

LargeDeformationJ2Plasticity::LargeDeformationJ2Plasticity(const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    _T(declareADProperty<Real>(prependBaseName("Temp"))),
    _T_old(getMaterialPropertyOldByName<Real>(prependBaseName("Temp"))),
    _sigma_0(getADMaterialProperty<Real>(prependBaseName("ref_yield_stress", true))),
    _sigma_y(declareADProperty<Real>(prependBaseName("yield_stress")))

{
}

void
LargeDeformationJ2Plasticity::initQpStatefulProperties()
{
  LargeDeformationPlasticityModel::initQpStatefulProperties();
  _sigma_y[_qp] = _sigma_0[_qp];
  _T[_qp] = 293;
}

void
LargeDeformationJ2Plasticity::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
{
  // First assume no plastic increment
  ADReal delta_ep = 0;
  _sigma_y[_qp] = _sigma_0[_qp] * (293 / _T[_qp]);

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

  // Update temp
  computeTemperature(delta_ep);
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
  return effective_trial_stress -
         _elasticity_model->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
             .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 1);
}

ADReal
LargeDeformationJ2Plasticity::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  return -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
              .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 2);
}

void
LargeDeformationJ2Plasticity::computeTemperature(const ADReal & delta_ep)
{

  _T[_qp] = (0.9 * _sigma_y[_qp] * delta_ep) / (8e-9 * 443e6) + _T_old[_qp];
}
