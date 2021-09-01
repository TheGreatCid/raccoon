//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "LargeDeformationPlasticityModel.h"

class LargeDeformationJ2Plasticity : public LargeDeformationPlasticityModel
{
public:
  static InputParameters validParams();

  LargeDeformationJ2Plasticity(const InputParameters & parameters);

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe) override;

protected:
  virtual void initQpStatefulProperties() override;

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ep) override;

  void computeTemperature(const ADReal & delta_ep, const ADReal & effective_stress);
  ADReal computeSigyDeriv(const ADReal & delta_ep,
                          const ADReal & effective_trial_stress,
                          const unsigned int derivative);
  MaterialProperty<Real> & _T;
  const MaterialProperty<Real> & _T_old;
  const MaterialProperty<Real> & _T0;
  const ADMaterialProperty<Real> & _sigma_0;
  const ADMaterialProperty<Real> & _n;
  const ADMaterialProperty<Real> & _ep0;
  const ADMaterialProperty<Real> & _xi;
  const ADMaterialProperty<Real> & _cv;
  const ADMaterialProperty<Real> & _rho;
  const ADMaterialProperty<Real> & _R;
};
