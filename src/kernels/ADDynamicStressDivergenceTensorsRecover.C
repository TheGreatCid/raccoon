//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDynamicStressDivergenceTensorsRecover.h"
#include "ADRankTwoTensorForward.h"
#include "ADReal.h"
#include "ElasticityTensorTools.h"

registerMooseObject("raccoonApp", ADDynamicStressDivergenceTensorsRecover);

InputParameters
ADDynamicStressDivergenceTensorsRecover::validParams()
{
  InputParameters params = ADStressDivergenceTensors::validParams();
  params.addClassDescription(
      "Residual due to stress related Rayleigh damping and HHT time integration terms.  "
      "When recompute_old_stress=true, sigma_old / sigma_older come from MOOSE's "
      "stateful-property tracked _stress_old / _stress_older everywhere -- which "
      "requires that ComputeLargeDeformationStress::initQpStatefulProperties is "
      "patched to evaluate the elasticity model on _Fm at INITIAL.  When false "
      "(default), sigma_old / sigma_older at t_step=1,2 are read from the legacy "
      "AD material properties `stress_sol` / `stress_old_store_sol`.");
  params.addParam<MaterialPropertyName>("zeta",
                                        0.0,
                                        "Name of material property or a constant real "
                                        "number defining the zeta parameter for the "
                                        "Rayleigh damping.");
  params.addParam<Real>("alpha", 0, "alpha parameter for HHT time integration");
  params.addParam<bool>("static_initialization",
                        false,
                        "Set to true to get the system to "
                        "equilibrium under gravity by running a "
                        "quasi-static analysis (by solving Ku = F) "
                        "in the first time step");
  params.addParam<bool>(
      "recompute_old_stress",
      false,
      "When true, sigma_old / sigma_older come from MOOSE's stateful-material "
      "_stress_old / _stress_older everywhere.  At t_step=1 _stress_old is the "
      "constitutive evaluated on F_recovered (via the patched "
      "ComputeLargeDeformationStress::initQpStatefulProperties), so sigma_old "
      "and sigma_curr share the same per-QP roundoff fingerprint and the HHT-alpha "
      "(1+alpha)*sigma - alpha*sigma_old cancellation is clean.  Requires the "
      "matching patches in ComputeLargeDeformationStress and ComputeDeformationGradient.  "
      "When false (default), reads sigma_old / sigma_older at t_step=1,2 from the "
      "legacy AD material properties `stress_sol` / `stress_old_store_sol`.");

  return params;
}

ADDynamicStressDivergenceTensorsRecover::ADDynamicStressDivergenceTensorsRecover(
    const InputParameters & parameters)
  : ADStressDivergenceTensors(parameters),
    _recompute_old_stress(getParam<bool>("recompute_old_stress")),
    _stress_older(getMaterialPropertyOlder<RankTwoTensor>(_base_name + "stress")),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>(_base_name + "stress")),
    _stress_older_sol(nullptr),
    _stress_old_sol(nullptr),
    _zeta(getMaterialProperty<Real>("zeta")),
    _alpha(getParam<Real>("alpha")),
    _static_initialization(getParam<bool>("static_initialization"))
{
  if (!_recompute_old_stress)
  {
    _stress_old_sol = &getADMaterialProperty<RankTwoTensor>("stress_sol");
    _stress_older_sol = &getADMaterialProperty<RankTwoTensor>("stress_old_store_sol");
  }
}

ADReal
ADDynamicStressDivergenceTensorsRecover::computeQpResidual()
{
  // Quasi-static first step: Div sigma only
  if (_static_initialization && _t == _dt)
  {
    ADReal residual = _stress[_qp].row(_component) * _grad_test[_i][_qp];
    if (_volumetric_locking_correction)
      residual +=
          _stress[_qp].trace() / 3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
    return residual;
  }

  if (_dt <= 0)
    return 0.0;

  // The old and older stress. Without recompute_old_stress, the stateful stress has no recovered
  // history on the first two steps, so the stress read from the recovery file stands in for it.
  ADRankTwoTensor stress_old = _stress_old[_qp];
  ADRankTwoTensor stress_older = _stress_older[_qp];
  if (!_recompute_old_stress && _t_step == 1)
  {
    stress_old = (*_stress_old_sol)[_qp];
    stress_older = (*_stress_older_sol)[_qp];
  }
  else if (!_recompute_old_stress && _t_step == 2)
    stress_older = (*_stress_old_sol)[_qp];

  // HHT-alpha + Rayleigh stress-divergence residual:
  //   R = [(1+alpha)*(1 + zeta/dt)]      * Div sigma
  //     - [alpha + (1+2*alpha)*zeta/dt]  * Div sigma_old
  //     + [alpha * zeta/dt]              * Div sigma_older
  const ADRankTwoTensor stress = (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) * _stress[_qp] -
                                 (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * stress_old +
                                 (_alpha * _zeta[_qp] / _dt) * stress_older;

  ADReal residual = stress.row(_component) * _grad_test[_i][_qp];
  if (_volumetric_locking_correction)
    residual += stress.trace() / 3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
  return residual;
}
