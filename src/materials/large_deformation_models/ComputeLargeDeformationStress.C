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

void
ComputeLargeDeformationStress::computeQpProperties()
{
  _elasticity_model->setQp(_qp);
  ADRankTwoTensor Fm_diff = _Fm[_qp] - (*_Fm_old)[_qp];

  double number_of_substeps = substepCheck(Fm_diff);
  if (number_of_substeps != 1)
  {
    std::cout << "get here" << std::endl;
    substepping(Fm_diff, number_of_substeps);
  }
  else
  {
    _elasticity_model->updateState(_Fm[_qp], _stress[_qp]);
  }
  if (_viscoelasticity_model)
  {
    _viscoelasticity_model->setQp(_qp);
    _stress[_qp] += _viscoelasticity_model->computeCauchyStress(_Fm[_qp], (*_Fm_old)[_qp]);
  }
}

double
ComputeLargeDeformationStress::substepCheck(ADRankTwoTensor & Fm_diff)
{
  // Will be defined as params later
  unsigned int number_of_substeps = 1;
  double max_inelastic_increment = 0.0001;
  //  unsigned int maximum_number_substeps = 25;
  double substep_strain_tolerance = 0.1;
  // Calculate Effective plastic strain to check for deformation tolerance
  const ADReal contracted_elastic_strain = (Fm_diff).doubleContraction(Fm_diff);
  const Real effective_elastic_strain =
      std::sqrt(3.0 / 2.0 * MetaPhysicL::raw_value(contracted_elastic_strain));
  const Real ratio = (effective_elastic_strain) / (max_inelastic_increment);
  if (ratio > substep_strain_tolerance)
    return number_of_substeps = std::ceil(ratio / substep_strain_tolerance);
  else
    return 1;
}

void
ComputeLargeDeformationStress::substepping(ADRankTwoTensor & Fm_diff, double & number_of_substeps)
{
  ADRankTwoTensor _temporary_deformation_gradient;
  Real dt_original = _dt;
  _dt = dt_original / number_of_substeps;
  std::cout << number_of_substeps << "--------------------------------" << std::endl;

  for (unsigned int istep = 0; istep < number_of_substeps; ++istep)
  {
    _temporary_deformation_gradient = (static_cast<Real>(istep) + 1) / number_of_substeps * Fm_diff;
    _temporary_deformation_gradient += _Fm[_qp];
    _elasticity_model->updateState(_temporary_deformation_gradient, _stress[_qp]);
  }
  // Reset Dt
  _dt = dt_original;
}
