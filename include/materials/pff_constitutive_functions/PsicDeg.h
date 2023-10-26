//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "PsicDegModel.h"

class PsicDeg : public PsicDegModel
{
public:
  static InputParameters validParams();

  PsicDeg(const InputParameters & parameters);

  virtual ADReal elasticEnergy(const ADReal & ep, const unsigned int derivative) override;
  virtual ADReal initialGuess(const ADReal & effective_trial_stress) override;

protected:
  const ADMaterialProperty<Real> & _psie_active;
  /// The fracture toughness
  const ADMaterialProperty<Real> & _Gc;

  /// The normalization constant
  const ADMaterialProperty<Real> & _c0;

  /// The regularization length
  const ADMaterialProperty<Real> & _l;

  const Real _eta;
  const Real _A;
  const Real _B;
  const Real _C;
  const Real _D;
  const VariableValue & _d;
  const Real _psic_orig;
};
