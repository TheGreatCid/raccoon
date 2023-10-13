//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADKernel.h"

class ADTimeDerivativeCustom : public ADKernel
{
public:
  static InputParameters validParams();

  ADTimeDerivativeCustom(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  const ADVariableValue & _d_dot;
};
