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

  /// Smooth Heaviside function for conditional logic
  ADReal smoothHeaviside(const ADReal x, const Real width) const;

  /// Smooth approximation of max(0, x)
  ADReal smoothMacaulay(const ADReal x, const Real eps) const;

  /// The penalty parameter
  const Real _penalty;

  /// Smoothing parameter
  const Real _epsilon;

  /// Previous time step value
  const VariableValue & _u_old;

  /// Use smooth approximation
  const bool _smooth;

  /// Use conditional bounds
  const bool _conditional;

  /// Enforce upper bound (true) or lower bound (false)
  const bool _upper;

  /// Upper bound value
  const Real _upper_val;

  /// Smoothing method
  const MooseEnum _smoothing_type;

  /// Adaptive penalty parameters
  const bool _adaptive_penalty;
  const Real _initial_penalty;
  const Real _penalty_growth;
  const Real _violation_tol;

  /// Smooth transition width for conditional mode
  const Real _transition_width;
};
