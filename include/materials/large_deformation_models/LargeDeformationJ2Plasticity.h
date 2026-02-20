//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "LargeDeformationPlasticityModel.h"
#include "RankTwoTensorForward.h"

class LargeDeformationJ2Plasticity : public LargeDeformationPlasticityModel
{
public:
  static InputParameters validParams();

  LargeDeformationJ2Plasticity(const InputParameters & parameters);

  virtual void updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe) override;

protected:
  virtual ADReal initialGuess(const ADReal & effective_trial_stress) override;

  virtual ADReal minimumPermissibleValue(const ADReal &) const override { return 0; }

  virtual ADReal computeResidual(const ADReal & effective_trial_stress,
                                 const ADReal & delta_ep) override;

  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & delta_ep) override;

  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & delta_ep) override;

  ADMaterialProperty<Real> & _phi;
  ADMaterialProperty<Real> & _flowstress;
  ADMaterialProperty<Real> & _visflowstress;

  /// When true, solve for (delta_ep, delta_T) together in a local 2x2 Newton
  const bool _coupled_temp_solve;
  /// Mass density material property (used when coupled_temp_solve = true)
  const ADMaterialProperty<Real> * _rho;
  /// Specific heat capacity c_v [J/(kg·K)] (used when coupled_temp_solve = true)
  const Real _cv;

  /// Converged dep from the 2x2 Newton; used to seed returnMappingSolve with a good initial guess
  Real _dep_2x2_guess;
  bool _use_2x2_guess;
};
