// JohnsonCookHardening.h
// This file is part of the RACCOON application
// being developed at Dolbow lab at Duke University
// http://dolbow.pratt.duke.edu

#pragma once

#include "PlasticHardeningModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"
#include <vector>

class JohnsonCookHardening : public PlasticHardeningModel,
                             public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  JohnsonCookHardening(const InputParameters & parameters);

  virtual ADReal initialGuess(const ADReal & effective_trial_stress) override;
  virtual ADReal plasticEnergy(const ADReal & ep, const unsigned int derivative = 0) override;
  virtual ADReal plasticDissipation(const ADReal & delta_ep,
                                    const ADReal & ep,
                                    const unsigned int derivative) override;
  virtual ADReal thermalConjugate(const ADReal & ep) override;

  // Set a per-QP local temperature (used inside a local Newton for this qp)
  virtual void setLocalTemperature(Real T) override;
  virtual void clearLocalTemperature() override;

  // Return the OLD timestep temperature at this QP (previous step)
  virtual Real getQpTemperatureOld() const override;

  // Aux helpers used by the plasticity material
  virtual Real temperatureDependenceLogDerivative(Real T) override;
  virtual Real thermalConjugateTemperatureDerivative(Real ep) override;
  virtual Real thermalConjugateEpDerivative(Real ep) override;
  virtual Real dissipationFlowStressRateJacobian(Real dep, Real ep) override;
  virtual ADReal getQpTAD() const override;

protected:
  // @{ The plastic energy parameters
  const ADMaterialProperty<Real> & _sigma_0;
  const Real _n;
  const Real _m;
  const Real _ep0;
  const Real _epdot0;
  const Real _T0;
  // Current iterative temperature (array per-qp) — AD for global Jacobian coupling
  const ADVariableValue & _T;
  // Old timestep temperature (array per-qp)
  const VariableValue & _T_old;
  const Real _tqf;
  const ADMaterialProperty<Real> & _A;
  const Real _B;
  const Real _C;
  const Real _Tm;
  // @}

  /// Name of the phase-field variable
  const VariableName _d_name;

  // @{ Plastic energy density and its derivative w/r/t damage
  const MaterialPropertyName _psip_name;
  ADMaterialProperty<Real> & _psip;
  ADMaterialProperty<Real> & _psip_active;
  ADMaterialProperty<Real> & _dpsip_dd;
  // @}

  // @{ The degradation function and its derivative w/r/t damage
  const MaterialPropertyName _gp_name;
  const ADMaterialProperty<Real> & _gp;
  const ADMaterialProperty<Real> & _dgp_dd;
  // @}

  /// Option to remove dissipation contributions
  const bool _disable_dissipation;

  /// Per-QP local temperature override for coupled thermo-plastic Newton solve
  /// Indexed by _qp in per-QP calls to setLocalTemperature/clearLocalTemperature/getQpT.
  std::vector<Real> _T_local;
  std::vector<char> _use_local_T;

private:
  ADReal temperatureDependence();

  // Return the temperature to use in material calls: local override (offset onto global AD T)
  // or the current global T. Returns ADReal so all material properties carry ∂/∂T_DOF.
  ADReal getQpT() const;
};
