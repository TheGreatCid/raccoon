//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputeLargeDeformationStress.h"
#include "LargeDeformationElasticityModel.h"
#include "LargeDeformationPlasticityModel.h"
#include "LargeDeformationViscoelasticityModel.h"

registerMooseObject("raccoonApp", ComputeLargeDeformationStress);

InputParameters
ComputeLargeDeformationStress::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("Stress calculator given an elasticity model, a plasticity model and "
                             "a viscoelasticity model. Large deformation is assumed.");

  params.addRequiredParam<MaterialName>("elasticity_model",
                                        "Name of the elastic stress-strain constitutive model");
  params.addParam<MaterialName>("plasticity_model", "Name of the plasticity model");
  params.addParam<MaterialName>("viscoelasticity_model", "Name of the viscoelasticity model");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

ComputeLargeDeformationStress::ComputeLargeDeformationStress(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _Fm(getADMaterialProperty<RankTwoTensor>(prependBaseName("mechanical_deformation_gradient"))),
    _Fm_old(isParamValid("viscoelasticity_model")
                ? &getMaterialPropertyOld<RankTwoTensor>(
                      prependBaseName("mechanical_deformation_gradient"))
                : nullptr),
    _stress(declareADProperty<RankTwoTensor>(prependBaseName("stress")))
{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}

void
ComputeLargeDeformationStress::initialSetup()
{
  _elasticity_model =
      dynamic_cast<LargeDeformationElasticityModel *>(&getMaterial("elasticity_model"));
  if (!_elasticity_model)
    paramError("elasticity_model",
               "Elasticity model " + _elasticity_model->name() +
                   " is not compatible with ComputeLargeDeformationStress");

  _plasticity_model =
      isParamValid("plasticity_model")
          ? dynamic_cast<LargeDeformationPlasticityModel *>(&getMaterial("plasticity_model"))
          : nullptr;
  if (_plasticity_model)
    _elasticity_model->setPlasticityModel(_plasticity_model);

  _viscoelasticity_model = isParamValid("viscoelasticity_model")
                               ? dynamic_cast<LargeDeformationViscoelasticityModel *>(
                                     &getMaterial("viscoelasticity_model"))
                               : nullptr;
}

void
ComputeLargeDeformationStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
}
// Total formulation for deformation gradient
void
ComputeLargeDeformationStress::computeQpProperties()
{
  _elasticity_model->setQp(_qp);

  // Call method to check for substepping
  bool substep = _elasticity_model->substepCheck(_Fm[_qp]);

  // Run substep if condition is true
  if (substep)
    substepping();

  _elasticity_model->updateState(_Fm[_qp], _stress[_qp]);

  if (_viscoelasticity_model)
  {
    _viscoelasticity_model->setQp(_qp);
    _stress[_qp] += _viscoelasticity_model->computeCauchyStress(_Fm[_qp], (*_Fm_old)[_qp]);
  }
}

// Substepping
void
ComputeLargeDeformationStress::substepping()
{
  //  std::cout << "substepping" << std::endl;
  // Store orignal _dt; Reset at the end of solve
  Real dt_original = _dt;
  // cut the original timestep
  _dt = dt_original / 100; // total_number_substeps;

  // Need algorithm to calculate this number
  Real numsubstep = 100;
  // initialize the inputs

  // substepping loop
  _elasticity_model->substepping(numsubstep, _Fm[_qp], _stress[_qp]);
  // Recover original timestep
  _dt = dt_original;
}
