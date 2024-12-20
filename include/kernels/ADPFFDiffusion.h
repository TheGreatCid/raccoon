//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADKernel.h"
#include "BaseNameInterface.h"
#include "libmesh/numeric_vector.h"

class ADPFFDiffusion : public ADKernel, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ADPFFDiffusion(const InputParameters & parameters);

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

  const VariableGradient & _d_old_grad;

  const VariableValue & _d_old_grad_ref;

  // MooseWritableVariable * _d_diff;
  // MooseVariable & _d_diff;

  // NumericVector<Real> & _d_diff;

  const VariableValue & _grad_dx;
  const VariableValue & _grad_dy;
  const VariableValue & _grad_dz;

  const VariableValue & _Fnobar_xx_1;
  const VariableValue & _Fnobar_xx_2;

  const VariableValue & _F_current;
};
