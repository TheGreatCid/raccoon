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
      "Selects how sigma_old / sigma_older are obtained at the first restart step via "
      "the `recompute_old_stress` flag (default true = constitutive-recomputed via MOOSE's "
      "_stress_old; false = read from `stress_sol` / `stress_old_store_sol` AD material "
      "properties as in the original SolutionTensor pattern).");
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
      true,
      "When true (default), sigma_old/sigma_older come from MOOSE's stateful-material "
      "tracked _stress_old / _stress_older.  At t_step=1 those are populated by the "
      "constitutive's INITIAL computeQpProperties() evaluation on F = F_recovered, so "
      "the per-QP roundoff fingerprint matches sigma_curr and the HHT-alpha "
      "(1+alpha)*sigma_curr - alpha*sigma_old cancellation is clean.  When false, the "
      "legacy SolutionTensor-fed path is used and `stress_sol` / `stress_old_store_sol` "
      "AD material properties must be declared in the input file.");
  // `solution` was required by the legacy implementation but never used in the
  // residual itself; keep it optional so recompute-mode input files don't need it.
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
    // Legacy mode: pull old/older sigma from SolutionTensor-style AD material
    // properties.  These must be declared in the input file.
    _stress_old_sol = &getADMaterialProperty<RankTwoTensor>("stress_sol");
    _stress_older_sol = &getADMaterialProperty<RankTwoTensor>("stress_old_store_sol");
  }

  // SolutionUserObject is optional in either mode (kernel doesn't use it).
  const auto solution_name = getParam<UserObjectName>("solution");
  if (!solution_name.empty())
    _solution_object_ptr = &getUserObject<SolutionUserObject>(solution_name);
}

ADReal
ADDynamicStressDivergenceTensorsRecover::computeQpResidual()
{
  /**
   * HHT-alpha + Rayleigh stress-divergence residual:
   *   R = [(1+alpha)*(1 + zeta/dt)] * Div sigma_curr
   *     - [alpha + (1+2*alpha)*zeta/dt] * Div sigma_old
   *     + [alpha * zeta/dt]            * Div sigma_older
   *
   * The branch structure handles two boundary conditions for the recover/
   * restart pattern: at t_step=1 (and t_step=2 in legacy mode) there is no
   * MOOSE-tracked older stress from a prior solve, so the source of
   * sigma_old / sigma_older is determined by `recompute_old_stress`.
   */
  ADReal residual;

  if (_static_initialization && _t == _dt)
  {
    // If static initialization is true, then in the first step residual is only Ku which is
    // stress.grad(test).
    residual = _stress[_qp].row(_component) * _grad_test[_i][_qp];

    if (_volumetric_locking_correction)
      residual +=
          _stress[_qp].trace() / 3.0 * (_avg_grad_test[_i] - _grad_test[_i][_qp](_component));
  }
  else if (_recompute_old_stress)
  {
    // RECOMPUTE PATH (default).
    //
    // At t_step=1 MOOSE's stateful-material machinery has populated _stress_old
    // from the INITIAL computeQpProperties() call of the constitutive, which
    // evaluated on F = F_recovered (seeded by ComputeDeformationGradient with
    // recover=true).  So _stress_old at t_step=1 IS the constitutive-recomputed
    // sigma at the recovered state.  _stress_older at t_step=1 is zero (no
    // prior state); the alpha*zeta/dt coefficient on it is the only place that
    // missing-history shows up, and that vanishes for zeta=0.
    //
    // No t_step branching is needed -- the same formula works at every step.
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
    // LEGACY PATH (recompute_old_stress = false).
    // First restart step: pull sigma_old from the SolutionTensor (= dumped sigma values).
    // sigma_older is the identity placeholder which is multiplied by alpha*zeta/dt and so
    // only contributes when Rayleigh damping is on.
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
  else if (_dt > 0 && _t_step == 2) // Need to use stored stress as older stress
  {
    // LEGACY PATH t_step=2: sigma_old is MOOSE-tracked (= the converged sigma at t_step=1);
    // sigma_older still comes from the SolutionTensor as the previous-previous slot.
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
    // LEGACY PATH t_step >= 3: fully MOOSE-tracked history.
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
