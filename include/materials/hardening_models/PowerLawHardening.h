//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "PlasticHardeningModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class PowerLawHardening : public PlasticHardeningModel,
                          public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  PowerLawHardening(const InputParameters & parameters);

  virtual ADReal plasticEnergy(const ADReal & ep, const unsigned int derivative = 0) override;
  virtual ADReal plasticDissipation(const ADReal & delta_ep,
                                    const ADReal & ep,
                                    const unsigned int derivative) override;
  ADReal piecewise();
  Real val(Real Tl, const std::vector<Real> & coeff, Real breaks);

protected:
  // @{ The plastic energy parameters
  const ADMaterialProperty<Real> & _n;
  const ADMaterialProperty<Real> & _ep0;
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
  const ADMaterialProperty<Real> & _sigma_0;
  ADMaterialProperty<Real> & _sigma_y;
  const Real _T0;
  const enum class Sigy_func { piece, exp, tan } _sigy_func;
  const ADVariableValue & _T;
  const Real _tqf;
  const Real _eps;

  // @}
};
