//* David devel

#pragma once

#include "LargeDeformationPlasticityModel.h"

class LargeDeformationJ2PlasticityModified : public LargeDeformationPlasticityModel
{
public:
  static InputParameters validParams();

  LargeDeformationJ2PlasticityModified(const InputParameters & parameters);

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe) override;

protected:
  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ep) override;
  virtual ADReal rexp_calc();

  virtual ADReal rtanh_calc();
  const ADVariableValue & _T;
  const Real _del;
  const Real _kappa;
  const Real _T0;
  const Real _ep0;
  const Real _ep0_dot;
  const Real _sigma0;
  const Real _k;
};
