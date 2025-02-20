//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "ADReal.h"
#include "EigenADReal.h"
#include "LargeDeformationJ2Plasticity.h"
#include "MooseError.h"
#include "RaccoonUtils.h"
#include <cmath>
#include <string>

registerMooseObject("raccoonApp", LargeDeformationJ2Plasticity);

InputParameters
LargeDeformationJ2Plasticity::validParams()
{
  InputParameters params = LargeDeformationPlasticityModel::validParams();
  params.addClassDescription("Large deformation $J_2$ plasticity. The exponential constitutive "
                             "update is used to update the plastic deformation.");
  params.addParam<bool>("recover", false, "do you want to recover");
  params.addParam<UserObjectName>("solution", "The SolutionUserObject to extract data from.");
  return params;
}

LargeDeformationJ2Plasticity::LargeDeformationJ2Plasticity(const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    _phi(declareADProperty<Real>("phi")),
    _flowstress(declareADProperty<Real>("flowstress")),
    _visflowstress(declareADProperty<Real>("visflowstress"))

{
  _check_range = true;

  if (!isParamValid("solution") && _recover == true)
    MaterialBase::mooseError("Need solution object!");

  if (_recover == true)
    _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");
}

void
LargeDeformationJ2Plasticity::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
{

  ADRankTwoTensor curr_Fp;

  // populate curr_FP and _ep_old_store
  //_ep_old store is a material property so that it can be used in the residual calculation

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

  if (_t_step == 4)
    std::cout << MetaPhysicL::raw_value(stress_dev_norm) << " " << delta_ep << std::endl;
  _phi[_qp] = computeResidual(stress_dev_norm, delta_ep);
  if (_phi[_qp] > 0)
    returnMappingSolve(stress_dev_norm, delta_ep, _console);
  if (_t_step == 4)
    std::cout << "after" << std::endl;
  _ep[_qp] = _ep_old[_qp] + delta_ep;

  if (_ep[_qp] == 0)
    _ep[_qp] = 1e-20;

  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  Fe = Fe * delta_Fp.inverse();
  stress = _elasticity_model->computeCauchyStress(Fe);

  _hardening_model->plasticEnergy(_ep[_qp]);
  _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 0);
  //  if(_t_step==1 && _ep[_qp]>0.001)
  //  	std::cout << _ep[_qp] << std::endl;

  if (_t_step > 0)
  {
    _heat[_qp] = _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 1) * delta_ep / _dt;

    _heat[_qp] += _hardening_model->thermalConjugate(_ep[_qp]) * delta_ep / _dt;
  }
  else
    _heat[_qp] = 0;

  _flowstress[_qp] = _hardening_model->plasticEnergy(_ep[_qp], 1);
  _visflowstress[_qp] = _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 1);
}

Real
LargeDeformationJ2Plasticity::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                       const ADReal & delta_ep)
{
  return raw_value(effective_trial_stress - _elasticity_model
                                                ->computeMandelStress(delta_ep * _Np[_qp],
                                                                      /*plasticity_update = */ true)
                                                .doubleContraction(_Np[_qp]));
}

ADReal
LargeDeformationJ2Plasticity::computeResidual(const ADReal & effective_trial_stress,
                                              const ADReal & delta_ep)
{
  ADReal ep;

  ep = _ep_old[_qp] + delta_ep;
  if (ep == 0)
  {
    ep = 1e-20;
  }
  if (_t_step == 4 && _qp == 1)
  {
    std::cout << "==========" << std::endl;
    std::cout << "delta_ep " << MetaPhysicL::raw_value(delta_ep) << std::endl;
    MetaPhysicL::raw_value(_Np[_qp]).print();

    // auto stress = _elasticity_model->computeMandelStress(delta_ep * _Np[_qp],
    //                                                      /*plasticity_update = */ true);

    // auto test = stress.doubleContraction(_Np[_qp]);
    // MetaPhysicL::raw_value(_Np[_qp]).print();
    // std::cout << "+++++++++++++++++++" << std::endl;
    // MetaPhysicL::raw_value(stress).print();
    // std::cout << "contraction " << MetaPhysicL::raw_value(test) << std::endl;

    std::cout << "==========" << std::endl;
  }
  // if (std::isnan(delta_ep))
  //   mooseError("ACHK");
  return effective_trial_stress -
         _elasticity_model
             ->computeMandelStress(delta_ep * _Np[_qp],
                                   /*plasticity_update = */ true)
             .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(ep, 1) -
         _hardening_model->plasticDissipation(delta_ep, ep, 1);
}

ADReal
LargeDeformationJ2Plasticity::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  ADReal ep;

  ep = _ep_old[_qp] + delta_ep;
  if (ep == 0)
  {
    ep = 1e-20;
  }
  return -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
              .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(ep, 2) -
         _hardening_model->plasticDissipation(delta_ep, ep, 2);
}
