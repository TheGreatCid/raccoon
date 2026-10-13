//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"

class SolutionUserObjectQP;

/**
 * A tensor material property read at each quadrature point from a recovery file
 */
class SolutionTensor : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  SolutionTensor(const InputParameters & parameters);

  void initStatefulProperties(unsigned int n_points) override;

  void initialSetup() override;

protected:
  const std::string _tensor_name;

  ADMaterialProperty<RankTwoTensor> & _tensor;
  /// Not read. Requesting the old value makes the property stateful, which is what gets
  /// initStatefulProperties() called and keeps the recovered value across time steps.
  const MaterialProperty<RankTwoTensor> & _tensor_old;
  const SolutionUserObjectQP * _solution_object_ptr;
};
