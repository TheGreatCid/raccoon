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
#include "Assembly.h"

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
  // `solution` is required by the legacy path but never used in the residual
  // itself; keep it optional so recompute-mode input files don't need it.
  params.addParam<UserObjectName>(
      "solution",
      "",
      "Optional SolutionUserObject.  Unused by the kernel itself; supported for "
      "input-file backwards compatibility with the legacy SolutionTensor pattern.");

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
    _static_initialization(getParam<bool>("static_initialization")),
    _solution_object_ptr(nullptr),
    _assembly_undisplaced(_fe_problem.assembly(_tid, this->_sys.number())),
    _q_point_undisplaced(_assembly_undisplaced.qPoints())
{
  if (!_recompute_old_stress)
  {
    _stress_old_sol = &getADMaterialProperty<RankTwoTensor>("stress_sol");
    _stress_older_sol = &getADMaterialProperty<RankTwoTensor>("stress_old_store_sol");
  }

  // getUserObject takes the PARAM NAME (not the resolved value) and looks up
  // the user object via that param.  Skip when the param value is empty so
  // recompute-mode inputs don't need to declare a SolutionUserObject.
  if (!getParam<UserObjectName>("solution").empty())
    _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");
}

ADReal
ADDynamicStressDivergenceTensorsRecover::computeQpResidual()
{
  /**
   * HHT-alpha + Rayleigh stress-divergence residual:
   *   R = [(1+alpha)*(1 + zeta/dt)] * Div sigma
   *     - [alpha + (1+2*alpha)*zeta/dt] * Div sigma_old
   *     + [alpha * zeta/dt]            * Div sigma_older
   */
  ADReal residual;

  if (_static_initialization && _t == _dt)
  {
    residual = _stress[_qp].row(_component) * _grad_test[_i][_qp];

    if (_volumetric_locking_correction)
      residual +=
          _stress[_qp].trace() / 3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
  }
  else if (_recompute_old_stress)
  {
    // RECOMPUTE PATH.  At t_step=1, _stress_old has been populated at INITIAL
    // by ComputeLargeDeformationStress::initQpStatefulProperties evaluating
    // the elasticity model on _Fm = F_recovered.  _stress_older at t_step=1
    // is zero (no prior state); its only coefficient is alpha*zeta/dt, which
    // vanishes for the zeta=0 (no Rayleigh) case.
    if (_dt > 0)
    {
      residual =
          _stress[_qp].row(_component) * _grad_test[_i][_qp] *
              (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
          (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * _stress_old[_qp].row(_component) *
              _grad_test[_i][_qp] +
          (_alpha * _zeta[_qp] / _dt) * _stress_older[_qp].row(_component) * _grad_test[_i][_qp];

      if (_volumetric_locking_correction)
        residual += (_stress[_qp].trace() * (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
                     (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * _stress_old[_qp].trace() +
                     (_alpha * _zeta[_qp] / _dt) * _stress_older[_qp].trace()) /
                    3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
    }
    else
      residual = 0.0;
  }
  else if (_dt > 0 && _t_step == 1)
  {
    // LEGACY PATH t_step=1: sigma_old from SolutionTensor (dumped sigma);
    // sigma_older is the identity placeholder, only contributes when zeta != 0.
    residual =
        _stress[_qp].row(_component) * _grad_test[_i][_qp] *
            (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
        (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * (*_stress_old_sol)[_qp].row(_component) *
            _grad_test[_i][_qp] +
        (_alpha * _zeta[_qp] / _dt) * (*_stress_older_sol)[_qp].row(_component) * _grad_test[_i][_qp];

    if (_volumetric_locking_correction)
      residual +=
          (_stress[_qp].trace() * (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
           (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * (*_stress_old_sol)[_qp].trace() +
           (_alpha * _zeta[_qp] / _dt) * (*_stress_older_sol)[_qp].trace()) /
          3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
  }
  else if (_dt > 0 && _t_step == 2)
  {
    // LEGACY PATH t_step=2: sigma_old MOOSE-tracked, sigma_older from SolutionTensor.
    residual =
        _stress[_qp].row(_component) * _grad_test[_i][_qp] *
            (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
        (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * _stress_old[_qp].row(_component) *
            _grad_test[_i][_qp] +
        (_alpha * _zeta[_qp] / _dt) * (*_stress_old_sol)[_qp].row(_component) * _grad_test[_i][_qp];

    if (_volumetric_locking_correction)
      residual += (_stress[_qp].trace() * (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
                   (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * _stress_old[_qp].trace() +
                   (_alpha * _zeta[_qp] / _dt) * (*_stress_old_sol)[_qp].trace()) /
                  3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
  }
  else if (_dt > 0)
  {
    // LEGACY PATH t_step >= 3: fully MOOSE-tracked.
    residual =
        _stress[_qp].row(_component) * _grad_test[_i][_qp] *
            (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
        (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * _stress_old[_qp].row(_component) *
            _grad_test[_i][_qp] +
        (_alpha * _zeta[_qp] / _dt) * _stress_older[_qp].row(_component) * _grad_test[_i][_qp];

    if (_volumetric_locking_correction)
      residual += (_stress[_qp].trace() * (1.0 + _alpha + (1.0 + _alpha) * _zeta[_qp] / _dt) -
                   (_alpha + (1.0 + 2.0 * _alpha) * _zeta[_qp] / _dt) * _stress_old[_qp].trace() +
                   (_alpha * _zeta[_qp] / _dt) * _stress_older[_qp].trace()) /
                  3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
  }
  else
    residual = 0.0;

  return residual;
}
