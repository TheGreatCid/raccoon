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
  params.addParam<Real>("substep_strain_tolerance",
                        0.1,
                        "Maximum ratio of the initial elastic strain increment at start of the "
                        "return mapping solve to the maximum inelastic strain allowable in a "
                        "single substep. Reduce this value to increase the number of substeps");
  params.addParam<unsigned>("maximum_number_substeps",
                            1000,
                            "The maximum number of substeps allowed before cutting the time step.");
  params.addParam<Real>("max_inelastic_increment",
                        1e-4,
                        "The maximum inelastic strain increment allowed in a time step");
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
    _stress(declareADProperty<RankTwoTensor>(prependBaseName("stress"))),
    _maximum_number_substeps(getParam<unsigned>("maximum_number_substeps")),
    _max_inelastic_increment(getParam<Real>("max_inelastic_increment")),
    _substep_strain_tolerance(getParam<Real>("substep_strain_tolerance")),
    _substep_its(declareADProperty<Real>(prependBaseName("substep_its"))),
    _effective_fm_diff(declareADProperty<Real>(prependBaseName("effective_fm_diff"))),
    _ratio(declareADProperty<Real>(prependBaseName("ratio")))

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

  // Calculate deformation gradient increment
  const ADRankTwoTensor Fm_diff = _Fm[_qp] - (*_Fm_old)[_qp];

  // Check number of substeps
  unsigned int number_of_substeps = substepCheck(Fm_diff);
  if (number_of_substeps != 1)
  {
    //    std::cout << "get here" << std::endl;
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
  //  unsigned int maximum_number_substeps = 25;
  // Calculate Effective plastic strain to check for deformation tolerance
  const ADReal contracted_elastic_diff = (Fm_diff).doubleContraction(Fm_diff);
  // const ADReal contracted_elastic_diff = 0;
  // Query the hardening model so not hard coded to J2
  const Real effective_elastic_diff = std::sqrt(3.0 / 2.0 * raw_value(contracted_elastic_diff));

  if (!MooseUtils::absoluteFuzzyEqual(effective_elastic_diff, 0.0))
  {
    _effective_fm_diff[_qp] = effective_elastic_diff;
    const Real ratiotest = (effective_elastic_diff) / (_max_inelastic_increment);
    _ratio[_qp] = ratiotest;

    if (ratiotest > _substep_strain_tolerance)
    {
      _effective_fm_diff[_qp] = effective_elastic_diff;
      //  std::cout << "here2" << std::endl;
      number_of_substeps = std::ceil(ratiotest / _substep_strain_tolerance);
      //    std::cout << number_of_substeps << std::endl;
    }
  }
  if (number_of_substeps > _maximum_number_substeps)
    mooseException("The number of substeps computed exceeds the maximum_number_substeps. The "
                   "system time step will be cut.");
  return number_of_substeps;
}

void
ComputeLargeDeformationStress::substepping(const ADRankTwoTensor & Fm_diff,
                                           unsigned int & number_of_substeps)
{
  ADRankTwoTensor _temporary_deformation_gradient;

  // store original dt
  Real dt_original = _dt;

  // Split dt by number of substeps
  _dt = dt_original / number_of_substeps;
  // std::cout << number_of_substeps << "--------------------------------" << std::endl;

  for (unsigned int istep = 0; istep < number_of_substeps; ++istep)
  {
    // Define deformation increment for current substep
    _temporary_deformation_gradient = (static_cast<Real>(istep) + 1) / number_of_substeps * Fm_diff;
    // std::cout << "temp1--" << MetaPhysicL::raw_value(_temporary_deformation_gradient(1, 1))
    //           << std::endl;

    // Add increment to orignal deformation gradient
    _temporary_deformation_gradient += (*_Fm_old)[_qp];
    // std::cout << "tem2--" << MetaPhysicL::raw_value(_temporary_deformation_gradient(1, 1))
    //           << std::endl;
    // std::cout << "current--" << MetaPhysicL::raw_value(_Fm[_qp](1, 1)) << std::endl;

    // state update with current deformation gradient
    _substep_its[_qp] += 1;
    _elasticity_model->updateState(_temporary_deformation_gradient, _stress[_qp]);
  }
  // Reset Dt
  _dt = dt_original;
}
