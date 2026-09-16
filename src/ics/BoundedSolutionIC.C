#include "BoundedSolutionIC.h"

#include <limits>

registerMooseObject("raccoonApp", BoundedSolutionIC);

InputParameters
BoundedSolutionIC::validParams()
{
  InputParameters params = SolutionIC::validParams();
  params.addParam<Real>(
      "lower_bound", 0.0, "Recovered values below this are clamped up to it (default 0).");
  params.addParam<Real>("upper_bound",
                        std::numeric_limits<Real>::max(),
                        "Recovered values above this are clamped down to it (default: no upper "
                        "clamp).");
  params.addClassDescription(
      "Reads an initial condition from a SolutionUserObject like SolutionIC, then clamps the "
      "recovered value into [lower_bound, upper_bound]. Useful on restart when the remapper can "
      "return a phase-field value slightly outside its physical range (e.g. small negative d).");
  return params;
}

BoundedSolutionIC::BoundedSolutionIC(const InputParameters & parameters)
  : SolutionIC(parameters),
    _lower_bound(getParam<Real>("lower_bound")),
    _upper_bound(getParam<Real>("upper_bound"))
{
  if (_lower_bound > _upper_bound)
    paramError("upper_bound",
               "upper_bound (",
               _upper_bound,
               ") must be >= lower_bound (",
               _lower_bound,
               ").");
}

Real
BoundedSolutionIC::value(const Point & p)
{
  return std::min(_upper_bound, std::max(_lower_bound, SolutionIC::value(p)));
}
