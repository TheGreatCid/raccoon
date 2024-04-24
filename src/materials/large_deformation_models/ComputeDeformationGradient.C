//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "ComputeDeformationGradient.h"
#include "Material.h"
#include "MooseTypes.h"
#include "SolutionUserObject.h"

registerADMooseObject("raccoonApp", ComputeDeformationGradient);

InputParameters
ComputeDeformationGradient::validParams()
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
  params.addParam<MaterialPropertyName>("F_store", "F_store");
  params.addParam<bool>("recover", false, "Are you trying to recover");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<UserObjectName>("solution", "The SolutionUserObject to extract data from.");
  return params;
}

ComputeDeformationGradient::ComputeDeformationGradient(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _coord_sys(_assembly.coordSystem()),
    _ndisp(coupledComponents("displacements")),
    _disp(adCoupledValues("displacements")),
    _grad_disp(adCoupledGradients("displacements")),
    _volumetric_locking_correction(getParam<bool>("volumetric_locking_correction") &&
                                   !this->isBoundaryMaterial()),
    _current_elem_volume(_assembly.elemVolume()),
    _F(declareADProperty<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _F_NoFbar(declareADProperty<RankTwoTensor>(prependBaseName("dg_noFbar"))),
    _F_store_Fbar(declareADProperty<RankTwoTensor>(prependBaseName("deformation_gradient_Fbar"))),
    _F_store_Fbar_old(
        getMaterialPropertyOld<RankTwoTensor>(prependBaseName("deformation_gradient_Fbar"))),
    _Fm(declareADProperty<RankTwoTensor>(prependBaseName("mechanical_deformation_gradient"))),
    _Fg_names(prependBaseName(
        getParam<std::vector<MaterialPropertyName>>("eigen_deformation_gradient_names"))),
    _Fgs(_Fg_names.size()),
    _recover(getParam<bool>("recover")),
    _solution_object_ptr(NULL)

{
  for (unsigned int i = 0; i < _Fgs.size(); ++i)
    _Fgs[i] = &Material::getADMaterialProperty<RankTwoTensor>(_Fg_names[i]);

  if (MaterialBase::getParam<bool>("use_displaced_mesh"))
    MaterialBase::paramError("use_displaced_mesh",
                             "The strain calculator needs to run on the undisplaced mesh.");
}

void
ComputeDeformationGradient::initialSetup()
{
  if (!isParamValid("solution") && _recover == true)
    MaterialBase::mooseError("Need solution object!");

  displacementIntegrityCheck();

  if (_recover == true)
    _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");
  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp.push_back(&_ad_zero);
    _grad_disp.push_back(&_ad_grad_zero);
  }

  // Apply volume averaging to inputed deformation gradient
}

void
ComputeDeformationGradient::displacementIntegrityCheck()
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != MaterialBase::_mesh.dimension())
    MaterialBase::paramError(
        "displacements",
        "The number of variables supplied in 'displacements' must match the mesh dimension.");

  // Don't use F-bar in 1D
  if (_ndisp == 1 && _volumetric_locking_correction)
    MaterialBase::paramError("volumetric_locking_correction",
                             "has to be set to false for 1-D problems.");

  // Check for RZ
  if (getBlockCoordSystem() == Moose::COORD_RZ && _ndisp != 2)
    MaterialBase::paramError(
        "displacements",
        "There must be two displacement variables provided, one in r-direction another in "
        "z-direction");
}

void
ComputeDeformationGradient::initStatefulProperties(unsigned int n_points)
{
  if (_recover == true)
  {
    ADReal ave_F_det_init = 0;

    // Temp holder for current element def grad
    ADRankTwoTensor curr_F;
    // Get average
    for (_qp = 0; _qp < n_points; ++_qp)
    {

      Point curr_Point = _q_point[_qp];

      // Populate tensor from solution object
      for (int i_ind = 0; i_ind < 3; i_ind++)
        for (int j_ind = 0; j_ind < 3; j_ind++)
        {
          curr_F(i_ind, j_ind) = _solution_object_ptr->pointValue(
              1, curr_Point, "dg_noFbar_" + std::to_string(i_ind) + std::to_string(j_ind), nullptr);
        }

      _F_store_Fbar[_qp] = curr_F;
      ave_F_det_init += curr_F.det() * _JxW[_qp] * _coord[_qp];
    }
    // Get averaged initial deformation tensor
    ave_F_det_init /= _current_elem_volume;

    for (_qp = 0; _qp < n_points; ++_qp)
    {
      // Populate current degrad tensor
      Point curr_Point = _q_point[_qp];
      for (int i_ind = 0; i_ind < 3; i_ind++)
        for (int j_ind = 0; j_ind < 3; j_ind++)
        {
          curr_F(i_ind, j_ind) = _solution_object_ptr->pointValue(
              1, curr_Point, "dg_noFbar_" + std::to_string(i_ind) + std::to_string(j_ind), nullptr);
        }
      // Store value
      _F_store_Fbar[_qp] *= std::cbrt(ave_F_det_init / curr_F.det());
    }
  }
}

void
ComputeDeformationGradient::initQpStatefulProperties()
{
  _F[_qp].setToIdentity();
  _Fm[_qp].setToIdentity();
}

void
ComputeDeformationGradient::computeQpProperties()
{
  _F_store_Fbar[_qp] = _F_store_Fbar_old[_qp];
}

ADReal
ComputeDeformationGradient::computeQpOutOfPlaneGradDisp()
{
  if (!MooseUtils::absoluteFuzzyEqual(_q_point[_qp](0), 0.0))
    return (*_disp[0])[_qp] / _q_point[_qp](0);
  else
    return 0.0;
}

void
ComputeDeformationGradient::computeProperties()
{

  ADReal ave_F_det = 0;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    ADRankTwoTensor A = ADRankTwoTensor::initializeFromRows(
        (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    if (_coord_sys == Moose::COORD_RZ)
      A(2, 2) = computeQpOutOfPlaneGradDisp();
    _F[_qp] = A;
    _F[_qp].addIa(1.0);
    // Multiply in old deformation
    if (_recover == true)
    {
      // Volume average the stored F

      _F[_qp] *= _F_store_Fbar[_qp];
    }
    if (_volumetric_locking_correction)
      ave_F_det += _F[_qp].det() * _JxW[_qp] * _coord[_qp];
  }

  if (_volumetric_locking_correction)
    ave_F_det /= _current_elem_volume;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    // Store _F[_qp] before volume averaging
    _F_NoFbar[_qp] = _F[_qp];
    if (_volumetric_locking_correction)
      _F[_qp] *= std::cbrt(ave_F_det / _F[_qp].det());

    // Remove the eigen deformation gradient
    ADRankTwoTensor Fg(ADRankTwoTensor::initIdentity);
    for (auto Fgi : _Fgs)
      Fg *= (*Fgi)[_qp];
    _Fm[_qp] = Fg.inverse() * _F[_qp];
  }
}
