//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "LargeDeformationElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"
#include "CNHElasticityInterface.h"

class CNHIsotropicElasticity : public LargeDeformationElasticityModel,
                               public DerivativeMaterialPropertyNameInterface,
                               public CNHElasticityInterface
{
public:
  static InputParameters validParams();

  CNHIsotropicElasticity(const InputParameters & parameters);

  virtual ADRankTwoTensor computeMandelStress(const ADRankTwoTensor & Fe,
                                              const bool plasticity_update = false) override;

  const ADMaterialProperty<Real> & getK() const override { return _K; }
  const ADMaterialProperty<Real> & getG() const override { return _G; }
  const ADMaterialProperty<Real> & getDegradation() const override { return _g; }
  const ADMaterialProperty<Real> & getDegradationDerivative() const override { return _dg_dd; }

  ADMaterialProperty<Real> & getPsie() override { return _psie; }
  ADMaterialProperty<Real> & getPsieActive() override { return _psie_active; }
  ADMaterialProperty<Real> & getDpsieDD() override { return _dpsie_dd; }

  ADReal volumetricKirchhoffPressure(const ADReal & J) const override;
  ADReal volumetricEnergy(const ADReal & J) const override;
  bool supportsEnergySplit() const override { return true; }

protected:
  // @{ Decomposition methods
  virtual ADRankTwoTensor computeMandelStressNoDecomposition(const ADRankTwoTensor & Fe,
                                                             const bool plasticity_update);
  virtual ADRankTwoTensor computeMandelStressVolDevDecomposition(const ADRankTwoTensor & Fe,
                                                                 const bool plasticity_update);
  // @}

  /**
   * Optionally add an un-degraded volumetric barrier that resists element inversion in fully
   * damaged elements. The barrier contributes an un-degraded energy 0.5*kb*(ln J)^2 and stress
   * kb*(ln J)*I, which diverge as J -> 0, so a (nearly) stiffness-less damaged element still resists
   * collapse/fold. It is applied only where the phase field d exceeds the given threshold.
   */
  void applyInversionBarrier(ADRankTwoTensor & stress,
                             const ADReal & J,
                             const bool plasticity_update);

  /// The bulk modulus
  const ADMaterialProperty<Real> & _K;

  /// The shear modulus
  const ADMaterialProperty<Real> & _G;

  /// Name of the phase-field variable
  const VariableName _d_name;

  // @{ Strain energy density and its derivative w/r/t damage
  const MaterialPropertyName _psie_name;
  ADMaterialProperty<Real> & _psie;
  ADMaterialProperty<Real> & _psie_active;
  ADMaterialProperty<Real> & _dpsie_dd;
  // @}

  // @{ The degradation function and its derivative w/r/t damage
  const MaterialPropertyName _g_name;
  const ADMaterialProperty<Real> & _g;
  const ADMaterialProperty<Real> & _dg_dd;
  // @}

  /// Decomposittion types
  const enum class Decomposition { none, spectral, voldev } _decomposition;

  // @{ Un-degraded inversion barrier for fully damaged elements
  /// Whether to add the inversion barrier (default off)
  const bool _use_inversion_barrier;
  /// Barrier stiffness coefficient kb
  const Real _inversion_barrier_coef;
  /// Only apply the barrier where the phase field exceeds this threshold
  const Real _inversion_barrier_d_threshold;
  /// The phase field value (to gate the barrier on damaged elements)
  const ADVariableValue & _d;
  // @}
};
