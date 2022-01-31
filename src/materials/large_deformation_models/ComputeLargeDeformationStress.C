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
  const ADRankTwoTensor Fm_diff = _Fm[_qp] - (*_Fm_old)[_qp];

  unsigned int number_of_substeps = substepCheck(Fm_diff);
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

unsigned int
ComputeLargeDeformationStress::substepCheck(const ADRankTwoTensor & Fm_diff)
{
  // Will be defined as params later
  unsigned int number_of_substeps = 1;
  //  Real max_inelastic_increment = 0.001;
  //  unsigned int maximum_number_substeps = 25;
  Real substep_strain_tolerance = 0.1;
  // Calculate Effective plastic strain to check for deformation tolerance
  const ADReal contracted_elastic_diff = (Fm_diff).doubleContraction(Fm_diff);
  // const ADReal contracted_elastic_diff = 0;
  const Real effective_elastic_diff = std::sqrt(3.0 / 2.0 * raw_value(contracted_elastic_diff));

  // There is an issue with the division of effective_elastic_diff and max_inelastic_incrememnt
  // If they are not within a certain order of magnitude of each other it affects convergence.
  // I have no idea how this is...
  // Perhaps the number is too small?

  if (!MooseUtils::absoluteFuzzyEqual(effective_elastic_diff, 0.0))
  {
    // When ratio is a certain value it is leads to convergence issues
    // This seems to be at ratio > 1e-9??
    const Real ratiotest = (10000) * (effective_elastic_diff); // / (10 * max_inelastic_increment);
    if (ratiotest > substep_strain_tolerance)
    {
      std::cout << "here2" << std::endl;
      number_of_substeps = std::ceil(ratiotest / substep_strain_tolerance);
    }
  }
  return number_of_substeps;
}

void
ComputeLargeDeformationStress::substepping(const ADRankTwoTensor & Fm_diff,
                                           unsigned int & number_of_substeps)
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
