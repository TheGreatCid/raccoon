//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "IncrementalDeformationGradient.h"

registerADMooseObject("raccoonApp", IncrementalDeformationGradient);

InputParameters
IncrementalDeformationGradient::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription(
      "This class computes the deformation gradient. Eigen deformation gradients are extracted "
      "from the total deformation gradient. The F-bar approach can optionally be used to correct "
      "volumetric locking.");

  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigen_deformation_gradient_names", {}, "List of eigen deformation gradients to be applied");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

IncrementalDeformationGradient::IncrementalDeformationGradient(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _f(declareADProperty<RankTwoTensor>(prependBaseName("incremental_deformation_gradient"))),
    _F(getADMaterialPropertyByName<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _F_old(getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("deformation_gradient")))
{
  if (getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh", "The strain calculator needs to run on the undisplaced mesh.");
}

void
IncrementalDeformationGradient::initQpStatefulProperties()
{
  _f[_qp].setToIdentity();
}

void
IncrementalDeformationGradient::computeProperties()
{
  _f[_qp] = _F[_qp] * _F_old[_qp].inverse();
}
