//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"

/**
 * This class computes the deformation gradient
 */
class IncrementalDeformationGradient : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  IncrementalDeformationGradient(const InputParameters & parameters);

  void computeProperties() override;

protected:
  virtual void initQpStatefulProperties() override;

  /// The total deformation gradient
  ADMaterialProperty<RankTwoTensor> & _f;
  const ADMaterialProperty<RankTwoTensor> & _F;
  const MaterialProperty<RankTwoTensor> & _F_old;
};
