//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ADKernel.h"
#include "BaseNameInterface.h"

class ADPenaltyConstraint : public ADKernel, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ADPenaltyConstraint(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// The fracture toughness
  const Real _penalty;
};
