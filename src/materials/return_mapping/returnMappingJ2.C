#include "RaccoonUtils.h"
#include "returnMappingJ2.h"
#include "Moose.h"
#include "MooseEnum.h"
#include "MooseObject.h"
#include "ConsoleStreamInterface.h"
#include "Conversion.h"
#include "MathUtils.h"

#include "DualRealOps.h"

#include <limits>
#include <string>
#include <cmath>
#include <memory>

InputParameters
returnMappingJ2::validParams()
{
  InputParameters params = emptyInputParameters();

  params.addParam<Real>(
      "relative_tolerance", 1e-8, "Relative convergence tolerance for Newton iteration");
  params.addParam<Real>(
      "absolute_tolerance", 1e-11, "Absolute convergence tolerance for Newton iteration");
  params.addParam<Real>("acceptable_multiplier",
                        10,
                        "Factor applied to relative and absolute "
                        "tolerance for acceptable convergence if "
                        "iterations are no longer making progress");
  params.addParam<Real>("lower_bound", 0, "lower bound for return mapping");
  return params;
}

returnMappingJ2::returnMappingJ2(const InputParameters & parameters)
  : _relative_tolerance(parameters.get<Real>("relative_tolerance")),
    _absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
    _acceptable_multiplier(parameters.get<Real>("acceptable_multiplier")),
    _num_resids(30),
    _residual_history(_num_resids, std::numeric_limits<Real>::max()),
    _iteration(0),
    _initial_residual(0.0),
    _residual(0.0),
    _max_its(1000),
    _lower_bound(parameters.get<Real>("lower_bound"))
{
}

void
returnMappingJ2::returnMappingSolveJ2(const ADReal & effective_trial_stress, ADReal & scalar)
{
  if (!newtonIterations(effective_trial_stress, scalar))
  {
    // If guess from NR goes below zero then start bilinear iterations
    bilinearIterations(effective_trial_stress, scalar);
  }
}

bool
returnMappingJ2::newtonIterations(const ADReal & effective_trial_stress, ADReal & scalar)
{
  Real reference_residual = computeReferenceResidual(effective_trial_stress, scalar);
  ADReal scalar_old = scalar;
  ADReal scalar_increment = 0;
  _residual = computeResidual(effective_trial_stress, scalar);
  while (_iteration < _max_its && !converged(_residual, reference_residual) &&
         !convergedAcceptable(_iteration, reference_residual))
  {
    _residual = computeResidual(effective_trial_stress, scalar);
    scalar_increment = -_residual / computeDerivative(effective_trial_stress, scalar);
    scalar = scalar_old + scalar_increment;
    // Catch the lower limit
    if (scalar < _lower_bound)
    {
      scalar = scalar_old;
      // Stop newton iterations go to bo bifurication method
      return false;
    }
    scalar_old = scalar;
    // Catch bad values -- Need to replace with solve state methods
    if (std::isnan(_residual) || std::isinf(MetaPhysicL::raw_value(_residual)))
      std::cout << "NAN or INF residual encountered" << std::endl;

    if (_iteration == _max_its)
      std::cout << "REACHED MAX ITS" << std::endl;
  }
  return true;
}

void
returnMappingJ2::bilinearIterations(const ADReal & effective_trial_stress, ADReal & scalar)
{
  ADReal lower = _lower_bound;
  ADReal upper = scalar;
  _iteration = 0;
  ADReal midpoint;
  if (computeResidual(effective_trial_stress, lower) *
          computeResidual(effective_trial_stress, upper) >
      0)
  {
    mooseException("Bisection Failed, no meaningful root found");
  }

  while (_iteration < _max_its && std::abs(upper - lower) > _absolute_tolerance)
  {
    // Use scalar value before NR went below zero, find midpoint
    midpoint = (lower + upper) / 2;

    if (computeResidual(effective_trial_stress, midpoint) *
            computeResidual(effective_trial_stress, lower) <
        0)
    {
      upper = midpoint;
    }
    else
    {
      lower = midpoint;
    }
    // Catch bad values
    if (std::isnan(_residual) || std::isinf(MetaPhysicL::raw_value(_residual)))
      mooseException("NAN or INF residual encountered");

    if (_iteration == _max_its)
      mooseException("REACHED MAX ITS");
  }
  scalar = midpoint;
}

// This is pulled fropm ADSingleVariableReturnMapping
bool
returnMappingJ2::convergedAcceptable(const unsigned int it, const Real reference)
{
  // Require that we have at least done _num_resids evaluations before we allow for
  // acceptable convergence
  if (it < _num_resids)
    return false;

  // Check to see whether the residual has dropped by convergence_history_factor over
  // the last _num_resids iterations. If it has (which means it's still making progress),
  // don't consider it to be converged within the acceptable limits.
  const Real convergence_history_factor = 10.0;
  if (std::abs(_residual * convergence_history_factor) <
      std::abs(_residual_history[(it + 1) % _num_resids]))
    return false;

  // Now that it's determined that progress is not being made, treat it as converged if
  // we're within the acceptable convergence limits
  return converged(_residual / _acceptable_multiplier, reference);
}

// This is pulled fropm ADSingleVariableReturnMapping
bool
returnMappingJ2::converged(const ADReal & ad_residual, const Real reference)
{
  const Real residual = MetaPhysicL::raw_value(ad_residual);
  return (std::abs(residual) <= _absolute_tolerance ||
          std::abs(residual / reference) <= _relative_tolerance);
}
