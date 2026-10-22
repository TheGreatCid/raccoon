//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SolutionTensor.h"
#include "SolutionUserObjectQP.h"

registerADMooseObject("raccoonApp", SolutionTensor);

InputParameters
SolutionTensor::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.addClassDescription(
      "A tensor material property read at each quadrature point from a recovery file. The "
      "property is named after the tensor with the suffix '_sol'.");
  params.addRequiredParam<UserObjectName>("solution",
                                          "The SolutionUserObjectQP to extract data from.");
  params.addParam<std::string>("tensor_name", "stress", "Name of tensor");
  params.addParam<std::string>(
      "file_tensor_name",
      "Name of the tensor in the recovery file, when it differs from 'tensor_name' (which names "
      "the property)");
  return params;
}

SolutionTensor::SolutionTensor(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _tensor_name(getParam<std::string>("tensor_name")),
    _file_tensor_name(isParamValid("file_tensor_name") ? getParam<std::string>("file_tensor_name")
                                                       : _tensor_name),
    _tensor(declareADProperty<RankTwoTensor>(prependBaseName(_tensor_name + "_sol"))),
    _tensor_old(getMaterialPropertyOld<RankTwoTensor>(prependBaseName(_tensor_name + "_sol"))),
    _solution_object_ptr(nullptr)
{
}

void
SolutionTensor::initialSetup()
{
  _solution_object_ptr = &getUserObject<SolutionUserObjectQP>("solution");
}

void
SolutionTensor::initStatefulProperties(unsigned int n_points)
{
  for (_qp = 0; _qp < n_points; ++_qp)
    _tensor[_qp] = _solution_object_ptr->qpTensor(_t, *_current_elem, _file_tensor_name, _qp);
}
