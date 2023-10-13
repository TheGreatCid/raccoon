//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "RateDegradationFunction.h"

registerMooseObject("raccoonApp", RateDegradationFunction);

InputParameters
RateDegradationFunction::validParams()
{
  InputParameters params = DegradationFunctionBase::validParams();
  params.addClassDescription("defines the power degradation function $g(d) = exp(-ddot*p)$.");

  params.set<std::string>("function") = "exp(-ddot*p)";

  const std::vector<std::string> default_params = {"p"};
  params.set<std::vector<std::string>>("parameter_names") = default_params;
  return params;
}

RateDegradationFunction::RateDegradationFunction(const InputParameters & parameters)
  : DegradationFunctionBase(parameters)
{
}
