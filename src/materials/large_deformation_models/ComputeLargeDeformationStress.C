//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputeLargeDeformationStress.h"
#include "EigenADReal.h"
#include "LargeDeformationElasticityModel.h"
#include "LargeDeformationPlasticityModel.h"
#include "LargeDeformationViscoelasticityModel.h"
#include "MooseError.h"

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
    _stress(declareADProperty<RankTwoTensor>(prependBaseName("stress"))),
    _stress_old(getMaterialPropertyOld<RankTwoTensor>("stress")),
    _stress_old_store(declareADProperty<RankTwoTensor>("stress_old_store"))
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
               "Elasticity model " + getParam<MaterialName>("elasticity_model") +
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
  // Elastic-only path: populate _stress at INITIAL by running the elasticity
  // model on _Fm.  For the recover/restart pattern, ComputeDeformationGradient
  // seeds _Fm = Fg^-1 * F_recovered at INITIAL, so _stress comes out as
  // sigma_recovered.  When MOOSE then copies _stress -> _stress_old for
  // t_step=1, downstream kernels (e.g. ADDynamicStressDivergenceTensorsRecover
  // with recompute_old_stress=true) see the correct sigma_old without
  // needing a parallel SolutionTensor pipeline.
  //
  // For plastic or viscoelastic configurations the constitutive's return
  // mapping needs its own stateful internal variables (plastic_strain,
  // be_bar, viscoelastic history) which are uninitialized at INITIAL.
  // Calling updateState here would invoke the iterative return mapping on
  // garbage state and can fail to converge (or worse, silently corrupt
  // the plastic state seed).  Fall back to the original zero-out so plastic
  // recovery uses its existing SolutionTensor-fed sigma_old path unchanged.
  if (_plasticity_model || _viscoelasticity_model)
  {
    _stress[_qp].zero();
    return;
  }
  _elasticity_model->setQp(_qp);
  _elasticity_model->updateState(_Fm[_qp], _stress[_qp]);
}

void
ComputeLargeDeformationStress::computeQpProperties()
{

  _elasticity_model->setQp(_qp);
  _elasticity_model->updateState(_Fm[_qp], _stress[_qp]);
  _stress_old_store[_qp] = _stress_old[_qp];

  if (_viscoelasticity_model)
  {
    _viscoelasticity_model->setQp(_qp);
    _stress[_qp] += _viscoelasticity_model->computeCauchyStress(_Fm[_qp], (*_Fm_old)[_qp]);
  }
}
