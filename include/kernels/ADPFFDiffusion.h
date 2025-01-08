//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADKernel.h"
#include "Adaptivity.h"
#include "BaseNameInterface.h"
#include "MooseTypes.h"
#include "RankTwoTensorForward.h"
#include "libmesh/numeric_vector.h"
#include "SolutionUserObject.h"

class ADPFFDiffusion : public ADKernel, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ADPFFDiffusion(const InputParameters & parameters);

  void initialSetup() override;

protected:
  virtual ADReal computeQpResidual() override;

  /// The fracture toughness
  const ADMaterialProperty<Real> & _Gc;

  /// The normalization constant
  const ADMaterialProperty<Real> & _c0;

  /// The regularization length
  const ADMaterialProperty<Real> & _l;

  // is this recovering?
  const bool _recover;

  // std::vector<const ADVariableGradient *> _grad_disp;

  // const VariableValue & _grad_xx;
  // const VariableValue & _grad_xy;
  // const VariableValue & _grad_xz;
  // const VariableValue & _grad_yx;
  // const VariableValue & _grad_yy;
  // const VariableValue & _grad_yz;
  // const VariableValue & _grad_zx;
  // const VariableValue & _grad_zy;
  // const VariableValue & _grad_zz;

  const SolutionUserObject * _solution_object_ptr;
};
