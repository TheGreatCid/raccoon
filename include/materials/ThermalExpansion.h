//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"
#include "Function.h"

/**
 * This class computes the deformation gradient
 */
class ThermalExpansion : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ThermalExpansion(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  virtual void computeQpProperties() override;

  // The thermal eigen deformation gradient
  ADMaterialProperty<RankTwoTensor> & _Fg;
  const MaterialProperty<RankTwoTensor> & _Fg_old;

  // The thermal expansion coefficient
  const Real _alpha;

  // The current temperature
  const ADVariableValue & _T;
  const VariableValue & _T_old;

  // The reference temperature
  const VariableValue & _T0;
};
