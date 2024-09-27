//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "ADReal.h"
#include "EigenADReal.h"
#include "LargeDeformationJ2Plasticity.h"
#include "MooseError.h"
#include "RaccoonUtils.h"
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
    _visflowstress(declareADProperty<Real>("visflowstress")),
    _ep_old_store(declareADProperty<Real>("ep_old_store")),
    _kirchoff(declareADProperty<RankTwoTensor>("kirchoff")),
    _recover(getParam<bool>("recover")),
    _solution_object_ptr(NULL)
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
  if (_t_step < 2 && _recover == true)
  {
    std::vector<std::string> indices = {"x", "y", "z"};

    _ep_old_store[_qp] =
        _solution_object_ptr->directValue(_current_elem, "ep_" + std::to_string(_qp + 1));
    for (int i_ind = 0; i_ind < 3; i_ind++)
      for (int j_ind = 0; j_ind < 3; j_ind++)
      {
        //  std::cout << i_ind << " " << j_ind << std::endl;
        curr_Fp(i_ind, j_ind) = _solution_object_ptr->directValue(
            _current_elem, "Fp_" + indices[i_ind] + indices[j_ind] + "_" + std::to_string(_qp + 1));
      }
  }
  // First assume no plastic increment
  ADReal delta_ep = 0;
  if (_recover == true && _t_step < 2)
  {
    Fe = Fe * curr_Fp.inverse();
  }
  else
  {
    Fe = Fe * _Fp_old[_qp].inverse();
  }

  stress = _elasticity_model->computeMandelStress(Fe);
  _kirchoff[_qp] = stress;
  // Compute the flow direction following the Prandtl-Reuss flow rule.
  // We guard against zero denominator.
  ADRankTwoTensor stress_dev = stress.deviatoric();
  ADReal stress_dev_norm = stress_dev.doubleContraction(stress_dev);
  if (MooseUtils::absoluteFuzzyEqual(stress_dev_norm, 0))
    stress_dev_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
  stress_dev_norm = std::sqrt(1.5 * stress_dev_norm);
  _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;

  // Return mapping
  if (_recover == true && _t_step == 0)
  {
    _ep[_qp] = _ep_old_store[_qp];
    _Fp[_qp] = curr_Fp;
  }
  else
  {
    _phi[_qp] = computeResidual(stress_dev_norm, delta_ep);
    if (_phi[_qp] > 0)
      returnMappingSolve(stress_dev_norm, delta_ep, _console);
  }

  // Using stored ep on first time step if recovering
  if (_t_step == 1 && _recover == true)
    _ep[_qp] = _ep_old_store[_qp] + delta_ep;
  else
    _ep[_qp] = _ep_old[_qp] + delta_ep;

  if (_ep[_qp] == 0)
  {
    _ep[_qp] = 1e-20;
  }
  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);

  // Using stored Fp on first time step if recovering
  if (_t_step < 2 && _recover == true)
  {
    if (_t_step == 0)
      _Fp[_qp] = curr_Fp;

    if (_t_step == 1)
      _Fp[_qp] = delta_Fp * curr_Fp;
  }
  else
    _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  Fe = Fe * delta_Fp.inverse();
  stress = _elasticity_model->computeCauchyStress(Fe);

  _hardening_model->plasticEnergy(_ep[_qp]);
  _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 0);

  if (_t_step > 0)
  {
    _heat[_qp] = _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 1) * delta_ep / _dt;

    _heat[_qp] += _hardening_model->thermalConjugate(_ep[_qp]) * delta_ep / _dt;
  }
  else
  {
    _heat[_qp] = 0;
  }
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

  // Using stored Fp on first time step if recovering
  if (_t_step == 1 && _recover == true)
    ep = _ep_old_store[_qp] + delta_ep;
  else
    ep = _ep_old[_qp] + delta_ep;
  if (ep == 0)
  {
    ep = 1e-20;
  }
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

  // Using stored Fp on first time step if recovering
  if (_t_step == 1 && _recover == true)
    ep = _ep_old_store[_qp] + delta_ep;
  else
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
