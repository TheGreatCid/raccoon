#pragma once

#include "BaseNameInterface.h"
#include "InputParameters.h"
#include "MooseTypes.h"
#include "DualRealOps.h"
class ConsoleStream;

class returnMappingJ2
{
public:
  static InputParameters validParams();

  returnMappingJ2(const InputParameters & parameters);
  virtual ~returnMappingJ2() {}

protected:
  void returnMappingSolveJ2(const ADReal & effective_trial_stress, ADReal & scalar);
  virtual ADReal computeResidual(const ADReal & effective_trial_stress, const ADReal & scalar) = 0;
  virtual ADReal computeDerivative(const ADReal & effective_trial_stress,
                                   const ADReal & scalar) = 0;
  virtual Real computeReferenceResidual(const ADReal & effective_trial_stress,
                                        const ADReal & scalar) = 0;
  bool converged(const ADReal & residual, const Real reference);

private:
  bool newtonIterations(const ADReal & effective_trial_stress, ADReal & scalar);
  void bilinearIterations(const ADReal & effective_trial_stress, ADReal & scalar);

  bool convergedAcceptable(const unsigned int it, const Real reference);

  /// Relative convergence tolerance
  Real _relative_tolerance;

  /// Absolute convergence tolerance
  Real _absolute_tolerance;

  /// Multiplier applied to relative and absolute tolerances for acceptable convergence
  Real _acceptable_multiplier;

  /// Number of residuals to be stored in history
  const std::size_t _num_resids;

  /// History of residuals used to check whether progress is still being made on decreasing the residual
  std::vector<Real> _residual_history;

  unsigned int _iteration;

  ///@{ Residual values, kept as members to retain solver state for summary outputting
  ADReal _initial_residual;
  ADReal _residual;
  ///@}

  const unsigned int _max_its;
};
