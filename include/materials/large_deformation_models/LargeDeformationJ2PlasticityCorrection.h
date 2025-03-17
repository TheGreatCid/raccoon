//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "LargeDeformationPlasticityModel.h"

class LargeDeformationJ2PlasticityCorrection : public LargeDeformationPlasticityModel
{
public:
  static InputParameters validParams();

  LargeDeformationJ2PlasticityCorrection(const InputParameters & parameters);

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe) override;

protected:
  virtual ADReal initialGuess(const ADReal & effective_trial_stress) override
  {
    return _hardening_model->initialGuess(effective_trial_stress);
  }

  virtual void initQpStatefulProperties() override;

  virtual ADReal minimumPermissibleValue(const ADReal &) const override { return 0; }

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ep) override;

  void computeCorrectionTerm(const ADRankTwoTensor & devbebar);

  // Bebar state variable
  ADMaterialProperty<RankTwoTensor> & _bebar;
  const MaterialProperty<RankTwoTensor> & _bebar_old;
  // ADMaterialProperty<RankTwoTensor> & _f;
  const ADMaterialProperty<Real> & _G;
  const ADMaterialProperty<Real> & _K;
  const MaterialProperty<RankTwoTensor> & _F_old;
  const ADMaterialProperty<RankTwoTensor> & _F;
};
