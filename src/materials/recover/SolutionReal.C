//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SolutionReal.h"
#include "SolutionUserObjectQP.h"

registerMooseObject("raccoonApp", SolutionReal);

InputParameters
SolutionReal::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.addClassDescription(
      "A scalar material property read at each quadrature point from a recovery file. The "
      "property is named after the scalar with the suffix '_sol'.");
  params.addRequiredParam<UserObjectName>("solution",
                                          "The SolutionUserObjectQP to extract data from.");
  params.addParam<std::string>("mat_name", "stress", "Name of the scalar");
  return params;
}

SolutionReal::SolutionReal(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _mat_name(getParam<std::string>("mat_name")),
    _mat(declareADProperty<Real>(prependBaseName(_mat_name + "_sol"))),
    _mat_old(getMaterialPropertyOld<Real>(prependBaseName(_mat_name + "_sol"))),
    _solution_object_ptr(nullptr)
{
}

void
SolutionReal::initialSetup()
{
  _solution_object_ptr = &getUserObject<SolutionUserObjectQP>("solution");
}

void
SolutionReal::initStatefulProperties(unsigned int n_points)
{
  for (_qp = 0; _qp < n_points; ++_qp)
    _mat[_qp] = _solution_object_ptr->qpValue(_t, *_current_elem, _mat_name, _qp);
}
