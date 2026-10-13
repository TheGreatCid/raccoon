//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "InputParameters.h"
#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"
#include "MaterialProperty.h"
#include "MooseArray.h"
#include "SolutionUserObject.h"
#include "Qp_Mapping.h"

/**
 * Per-QP scalar material that reads values from a SolutionUserObject (typical
 * use: the recovery dump).
 */
class SolutionReal : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  SolutionReal(const InputParameters & parameters);

  void initStatefulProperties(unsigned int n_points) override;

  void initialSetup() override;

protected:
  const std::string _mat_name;

  ADMaterialProperty<Real> & _mat;
  /// Not read. Requesting the old value makes the property stateful, which is what gets
  /// initStatefulProperties() called and keeps the recovered value across time steps.
  const MaterialProperty<Real> & _mat_old;

  const SolutionUserObject * _solution_object_ptr;

  QpMapping::Element _element;

  unsigned int _qpnum;

private:
  const std::unordered_map<int, int> * _lookup;
};
