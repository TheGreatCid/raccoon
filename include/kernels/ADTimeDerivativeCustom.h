//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADKernelValue.h"

class ADTimeDerivativeCustom : public ADKernelValue
{
public:
  static InputParameters validParams();

  ADTimeDerivativeCustom(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;
  const ADVariableValue _d_dot;
};
