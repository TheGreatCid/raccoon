//* David Development

#pragma once

#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"

class EffectiveStress : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  EffectiveStress(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  // @{ Decomposition methods
  // @}

  /// The bulk modulus
  const ADMaterialProperty<RankTwoTensor> & _stress;
  ADMaterialProperty<ADReal> & _effective_stress;
};
