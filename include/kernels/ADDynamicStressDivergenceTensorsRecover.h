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
 * Two paths for the old/older stress on the first restart step (when MOOSE's
 * own _stress_old would otherwise hold no usable value):
 *
 *  - recompute_old_stress = true  (default, recommended)
 *      Use MOOSE's _stress_old / _stress_older everywhere.  At t_step=1 the
 *      stateful-material machinery has populated _stress_old from the INITIAL
 *      computeQpProperties() call of the constitutive, which evaluates on
 *      F = F_recovered (seeded by ComputeDeformationGradient(recover=true)).
 *      σ_old therefore comes through the SAME constitutive code path as
 *      σ_curr, so their per-QP roundoff fingerprints match and the
 *      HHT-alpha (1+α)σ − α·σ_old cancellation is clean.
 *
 *  - recompute_old_stress = false  (legacy)
 *      Read σ_old / σ_older at t_step=1,2 from the AD material properties
 *      `stress_sol` / `stress_old_store_sol` -- typically `SolutionTensor`
 *      materials reading the dumped σ values from the reference's exodus
 *      file.  σ_old then carries the reference-side roundoff fingerprint;
 *      the kernel's HHT cancellation leaks K-amplified per-QP noise into
 *      the residual and the disp solution.
 */
class ADDynamicStressDivergenceTensorsRecover : public ADStressDivergenceTensors
{
public:
  static InputParameters validParams();

  ADDynamicStressDivergenceTensorsRecover(const InputParameters & parameters);

protected:
  ADReal computeQpResidual();

  /// If true, σ_old / σ_older come from MOOSE's stateful-property machinery
  /// (constitutive-recomputed); if false, from the legacy `stress_sol` /
  /// `stress_old_store_sol` AD material properties.
  const bool _recompute_old_stress;

  ///{@ MOOSE-tracked old/older stress (always bound; the only source when
  ///   _recompute_old_stress is true).
  const MaterialProperty<RankTwoTensor> & _stress_older;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  ///@}

  ///{@ Legacy SolutionTensor-fed old/older stress.  Only bound when
  ///   _recompute_old_stress == false.  Held as pointers so a recompute-mode
  ///   input file doesn't have to declare the `stress_sol` /
  ///   `stress_old_store_sol` materials.
  const ADMaterialProperty<RankTwoTensor> * _stress_older_sol;
  const ADMaterialProperty<RankTwoTensor> * _stress_old_sol;
  ///@}

  // Rayleigh damping parameter _zeta and HHT time integration parameter _alpha
  const MaterialProperty<Real> & _zeta;
  const Real _alpha;
  const bool _static_initialization;

  /// Optional SolutionUserObject -- never used in computeQpResidual, kept for
  /// backwards-compatible input-file syntax.  May be null in recompute mode.
  const SolutionUserObject * _solution_object_ptr;
  Assembly & _assembly_undisplaced;
  const MooseArray<Point> & _q_point_undisplaced;
};
