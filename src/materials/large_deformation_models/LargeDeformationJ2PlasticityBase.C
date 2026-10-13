//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "LargeDeformationJ2PlasticityBase.h"
#include "SolutionUserObjectQP.h"

InputParameters
LargeDeformationJ2PlasticityBase::validParams()
{
  return LargeDeformationPlasticityModel::validParams();
}

LargeDeformationJ2PlasticityBase::LargeDeformationJ2PlasticityBase(
    const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters)
{
}

Real
LargeDeformationJ2PlasticityBase::recoverNonNegative(const std::string & property) const
{
  return std::max(0., _solution_object_ptr->qpValue(_t, *_current_elem, property, _qp));
}

RankTwoTensor
LargeDeformationJ2PlasticityBase::recoverTensor(const std::string & property) const
{
  return _solution_object_ptr->qpTensor(_t, *_current_elem, property, _qp);
}

void
LargeDeformationJ2PlasticityBase::initQpStatefulProperties()
{
  _Fp[_qp].setToIdentity();
  _ep[_qp] = 0;

  if (!_recover)
    return;

  _ep[_qp] = recoverNonNegative("effective_plastic_strain");
  _Fp[_qp] = recoverTensor("plastic_deformation_gradient");
}
