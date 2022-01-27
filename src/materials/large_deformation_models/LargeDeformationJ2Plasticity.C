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
  return params;
}

LargeDeformationJ2Plasticity::LargeDeformationJ2Plasticity(const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters), _heat(declareADProperty<Real>("heat"))
{
}
// Add method to check for substepping, call method from stress calculator
bool
LargeDeformationJ2Plasticity::substepCheck(ADRankTwoTensor & Fe)
{
  // Perform trial stress calculation
  ADReal delta_ep = 0;
  // Hold local form of stress
  ADRankTwoTensor stress = _elasticity_model->computeMandelStress(Fe);

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

  // Placeholder value. Need to calculate what is a reasonable threshhold
  if (phi > 1)
    return true;
  else
    return false;
}

ADReal
LargeDeformationJ2Plasticity::trialStress(ADRankTwoTensor & stress)
{
  ADRankTwoTensor stress_dev = stress.deviatoric();
  ADReal stress_dev_norm = stress_dev.doubleContraction(stress_dev);
  if (MooseUtils::absoluteFuzzyEqual(stress_dev_norm, 0))
    stress_dev_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
  stress_dev_norm = std::sqrt(1.5 * stress_dev_norm);
  _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;
  return stress_dev_norm;
}

void
LargeDeformationJ2Plasticity::substepping(ADReal numsubstep,
                                          ADRankTwoTensor & Fe,
                                          ADRankTwoTensor & /*stress*/)
{
  // Assumpe delta_ep = 0

  ADReal delta_ep = 0;
  ADRankTwoTensor stresslocal;
  Fe = Fe * _Fp_old[_qp].inverse();
  stresslocal = _elasticity_model->computeMandelStress(Fe);

  // Calculate stress where phi = 0
  ADReal effective_trial_stress = trialStress(stresslocal);
  ADReal phi = computeResidual(effective_trial_stress, delta_ep);

  ADReal base_trial_stress = phi;
  // Difference of trial stress from onset of plasticity to final stress
  ADReal stress_inc = base_trial_stress - effective_trial_stress;
  // Split the stress by number of substeps
  ADReal stress_step = stress_inc / numsubstep;

  // Begin substep look and feed split stress into return mapping
  for (unsigned int step = 0; step < 100 /*total_number_substeps*/; ++step)
  {
    // Increment stress by substep
    stresslocal = base_trial_stress + stress_step * step;
    effective_trial_stress = trialStress(stresslocal);

    ADReal phi = computeResidual(effective_trial_stress, delta_ep);

    if (phi > 0)
    {
      returnMappingSolve(effective_trial_stress, delta_ep, _console);
    }
    // Update plastic strain and heat
    _ep[_qp] = _ep_old[_qp] + delta_ep;

    _heat[_qp] = _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 1) * delta_ep / _dt;
    _heat[_qp] += _hardening_model->thermalConjugate(_ep[_qp]) * delta_ep / _dt;
    _hardening_model->plasticEnergy(_ep[_qp]);
  }
  // Updates deformation after return mapping and update heat
  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  // Update stress and energy
  updateState();
}
void
LargeDeformationJ2Plasticity::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
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
  {
    returnMappingSolve(stress_dev_norm, delta_ep, _console);
  }
  _ep[_qp] = _ep_old[_qp] + delta_ep;

  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  // Update stress and energy
  Fe = Fe * delta_Fp.inverse();
  stress = _elasticity_model->computeCauchyStress(Fe);
  _hardening_model->plasticEnergy(_ep[_qp]);

  // PowerLawHardening* ourHardeningCase = static_cast<PowerLawHardening*>(_hardening_model);

  // Compute generated heat
  // _heat[_qp] = ourHardeningCase->plasticDissipation(delta_ep, 1) * delta_ep / _dt;
  //  std::cout <<"0 - " <<raw_value(_heat[_qp]) << std::endl;

  _heat[_qp] = _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 1) * delta_ep / _dt;
  //  std::cout <<"1 - " <<raw_value(_heat[_qp]) << std::endl;
  _heat[_qp] += _hardening_model->thermalConjugate(_ep[_qp]) * delta_ep / _dt;
  // std::cout <<"2 - " <<raw_value(_heat[_qp]) << std::endl;
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
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 1) -
         _hardening_model->plasticDissipation(delta_ep, _ep_old[_qp] + delta_ep, 1);
}

ADReal
LargeDeformationJ2Plasticity::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  return -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
              .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(_ep_old[_qp] + delta_ep, 2) -
         _hardening_model->plasticDissipation(delta_ep, _ep_old[_qp] + delta_ep, 2);
}
