//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "PowerDegradationFunctionEtaMat.h"

registerMooseObject("raccoonApp", PowerDegradationFunctionEtaMat);

InputParameters
PowerDegradationFunctionEtaMat::validParams()
{
  InputParameters params = DegradationFunctionBase::validParams();
  params.addClassDescription(
      "defines the power degradation function $g(d) = (1-d)^p (1-\\eta) + \\eta$.");

  params.set<std::string>("expression") = "(1-d)^p*(1-eta)+eta";
  params.set<std::vector<std::string>>("parameter_names") = {"p"};
  params.set<std::vector<std::string>>("material_property_names") = {"eta"};

  return params;
}

PowerDegradationFunctionEtaMat::PowerDegradationFunctionEtaMat(const InputParameters & parameters)
  : DegradationFunctionBase(parameters)
{
}
