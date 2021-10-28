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
  params.addRequiredParam<MaterialPropertyName>("rho", "the density");
  params.addRequiredParam<MaterialPropertyName>("R", "TS_Rate");
  params.addRequiredParam<MaterialPropertyName>("c_v", "c_v");
  return params;
}

LargeDeformationJ2Plasticity::LargeDeformationJ2Plasticity(const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    _sigma_0(getADMaterialProperty<Real>(prependBaseName("ref_yield_stress"))),
    _n(getADMaterialProperty<Real>(prependBaseName("exponent", true))),
    _ep0(getADMaterialProperty<Real>(prependBaseName("reference_plastic_strain", true))),
    _xi(getADMaterialProperty<Real>(prependBaseName("Taylor_Quinney"))),
    _cv(getADMaterialProperty<Real>(prependBaseName("c_v"))),
    _rho(getADMaterialProperty<Real>(prependBaseName("rho"))),
    _R(getADMaterialProperty<Real>(prependBaseName("R"))),
    _delta_ep(declareADProperty<Real>("delta_ep")),
    _stress_eff(declareADProperty<Real>("stress_eff"))
{
}

void
LargeDeformationJ2Plasticity::initQpStatefulProperties()
{
  //_T[_qp] = _Tinit[_qp];
  LargeDeformationPlasticityModel::initQpStatefulProperties();
  //_sigma_y[_qp] = _sigma_0[_qp];
  // std::cout << _T0[_qp]<< std::endl;
}

void
LargeDeformationJ2Plasticity::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
{
  // First assume no plastic increment
  _delta_ep[_qp] = 0;
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
  ADReal phi = computeResidual(stress_dev_norm, _delta_ep[_qp]);
  if (phi > 0)
    returnMappingSolve(stress_dev_norm, _delta_ep[_qp], _console);
  _ep[_qp] = _ep_old[_qp] + _delta_ep[_qp];
  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(_delta_ep[_qp] * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  // Update stress
  Fe = Fe * delta_Fp.inverse();
  stress = _elasticity_model->computeCauchyStress(Fe);
  _hardening_model->plasticEnergy(_ep[_qp]);

  // Declare Effective Plastic Stress to output so can use in temp solver
  stress_dev = stress.deviatoric();
  stress_dev_norm = stress_dev.doubleContraction(stress_dev);
  _stress_eff[_qp] = std::sqrt(1.5 * stress_dev_norm);
}

Real
LargeDeformationJ2Plasticity::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                       const ADReal & _delta_ep)
{
  return MetaPhysicL::raw_value(effective_trial_stress -
                                _elasticity_model
                                    ->computeMandelStress(_delta_ep * _Np[_qp],
                                                          /*plasticity_update = */ true)
                                    .doubleContraction(_Np[_qp]));
}

ADReal
LargeDeformationJ2Plasticity::computeResidual(const ADReal & effective_trial_stress,
                                              const ADReal & _delta_ep)
{

  return effective_trial_stress -
         _elasticity_model->computeMandelStress(_delta_ep * _Np[_qp], /*plasticity_update = */ true)
             .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + _delta_ep, 1);
}

ADReal
LargeDeformationJ2Plasticity::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & _delta_ep)
{
  return -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
              .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + _delta_ep, 2);
}
