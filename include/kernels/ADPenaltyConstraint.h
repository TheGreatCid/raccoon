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
  const Real _epsilon;

  const VariableValue & _u_old;

  const bool _smooth;

  const bool _conditional;

  const bool _upper;

  const Real _upper_val;
};
