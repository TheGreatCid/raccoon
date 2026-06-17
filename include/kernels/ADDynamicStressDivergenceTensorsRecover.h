//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADStressDivergenceTensors.h"
#include "SolutionUserObject.h"

/**
 * AD HHT-alpha (+ Rayleigh) stress-divergence kernel for the recover/restart
 * pattern.
 *
 * Two paths for sigma_old / sigma_older on the first restart step, selected
 * via the `recompute_old_stress` flag:
 *
 *  - false (default, legacy):
 *      At t_step=1,2 read sigma_old / sigma_older from the AD material
 *      properties `stress_sol` / `stress_old_store_sol` -- typically
 *      SolutionTensor materials fed by a SolutionUserObject reading the
 *      reference dump's exodus values.
 *
 *  - true (requires the patched
 *      ComputeLargeDeformationStress::initQpStatefulProperties that runs
 *      the elasticity model on _Fm at INITIAL, plus the matching
 *      ComputeDeformationGradient::initStatefulProperties patch that seeds
 *      _Fm = Fg^-1 * F_recovered):
 *      Use MOOSE's stateful-material _stress_old / _stress_older
 *      everywhere.  At t_step=1, _stress_old is the constitutive evaluated
 *      on F_recovered -- so sigma_old comes through the same code path as
 *      sigma_curr and the per-QP roundoff fingerprints match.  No
 *      SolutionTensor materials are needed.
 */
class ADDynamicStressDivergenceTensorsRecover : public ADStressDivergenceTensors
{
public:
  static InputParameters validParams();

  ADDynamicStressDivergenceTensorsRecover(const InputParameters & parameters);

protected:
  ADReal computeQpResidual();

  /// If true, use MOOSE-tracked _stress_old / _stress_older everywhere.
  const bool _recompute_old_stress;

  ///{@ MOOSE-tracked old/older stress (always bound).
  const MaterialProperty<RankTwoTensor> & _stress_older;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  ///@}

  ///{@ Legacy SolutionTensor-fed old/older stress.  Only bound when
  ///   _recompute_old_stress == false so recompute-mode input files don't
  ///   need to declare `stress_sol` / `stress_old_store_sol`.
  const ADMaterialProperty<RankTwoTensor> * _stress_older_sol;
  const ADMaterialProperty<RankTwoTensor> * _stress_old_sol;
  ///@}

  // Rayleigh damping parameter _zeta and HHT time integration parameter _alpha
  const MaterialProperty<Real> & _zeta;
  const Real _alpha;
  const bool _static_initialization;

  /// Optional SolutionUserObject -- never used in computeQpResidual.  May
  /// be null in recompute mode.
  const SolutionUserObject * _solution_object_ptr;
  Assembly & _assembly_undisplaced;
  const MooseArray<Point> & _q_point_undisplaced;
};
