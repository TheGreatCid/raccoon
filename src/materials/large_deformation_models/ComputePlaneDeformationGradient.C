//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputePlaneDeformationGradient.h"
#include "Qp_Mapping.h"

registerADMooseObject("raccoonApp", ComputePlaneDeformationGradient);

InputParameters
ComputePlaneDeformationGradient::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription(
      "This class computes the deformation gradient. Eigen deformation gradients are extracted "
      "from the total deformation gradient. The F-bar approach can optionally be used to correct "
      "volumetric locking.");
  params.addCoupledVar("out_of_plane_strain", "strain_zz");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigen_deformation_gradient_names", {}, "List of eigen deformation gradients to be applied");
  params.addParam<bool>("recover", false, "Are you trying to recover");
  params.addParam<UserObjectName>("solution", "The SolutionUserObject to extract data from.");
  params.addParam<Real>("num_qps", 8, "Number of QPs");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addRequiredParam<MooseEnum>(
      "element", MooseEnum(QpMapping::ELEMENT_ENUM_DEFINITION), "The element type");
  return params;
}

ComputePlaneDeformationGradient::ComputePlaneDeformationGradient(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _coord_sys(_assembly.coordSystem()),
    _out_of_plane_strain(adCoupledValue("out_of_plane_strain")),
    _ndisp(coupledComponents("displacements")),
    _disp(adCoupledValues("displacements")),
    _grad_disp(adCoupledGradients("displacements")),
    _volumetric_locking_correction(getParam<bool>("volumetric_locking_correction") &&
                                   !this->isBoundaryMaterial()),
    _current_elem_volume(_assembly.elemVolume()),
    _F(declareADProperty<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _Fnobar(declareADProperty<RankTwoTensor>(prependBaseName("Fnobar"))),
    _F_store_Fbar(declareADProperty<RankTwoTensor>(prependBaseName("deformation_gradient_Fbar"))),
    _F_store_Fbar_old(
        getMaterialPropertyOld<RankTwoTensor>(prependBaseName("deformation_gradient_Fbar"))),
    _Fm(declareADProperty<RankTwoTensor>(prependBaseName("mechanical_deformation_gradient"))),
    _F_store_noFbar(
        declareADProperty<RankTwoTensor>(prependBaseName("deformation_gradient_NoFbar"))),
    _F_store_noFbar_old(
        getMaterialPropertyOld<RankTwoTensor>(prependBaseName("deformation_gradient_NoFbar"))),
    _Fg_names(prependBaseName(
        getParam<std::vector<MaterialPropertyName>>("eigen_deformation_gradient_names"))),
    _Fgs(_Fg_names.size()),
    _recover(getParam<bool>("recover")),
    _solution_object_ptr(NULL),
    _element(getParam<MooseEnum>("element").getEnum<QpMapping::Element>()),
    _Frobenius(declareProperty<Real>(prependBaseName("Frobenius_norm")))

{
  for (unsigned int i = 0; i < _Fgs.size(); ++i)
    _Fgs[i] = &getADMaterialProperty<RankTwoTensor>(_Fg_names[i]);

  if (getParam<bool>("use_displaced_mesh"))
    paramError("use_displaced_mesh", "The strain calculator needs to run on the undisplaced mesh.");
  if (_recover)
    _lookup = QpMapping::getLookup(_element, _qpnum, /*reversed=*/true);
}

void
ComputePlaneDeformationGradient::initialSetup()
{
  displacementIntegrityCheck();

  // set unused dimensions to zero
  for (unsigned i = _ndisp; i < 3; ++i)
  {
    _disp.push_back(&_ad_zero);
    _grad_disp.push_back(&_ad_grad_zero);
  }
  if (_recover == true)
    _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");
}

void
ComputePlaneDeformationGradient::initStatefulProperties(unsigned int n_points)
{
  using std::cbrt;

  for (_qp = 0; _qp < n_points; ++_qp)
  {
    _F[_qp].setToIdentity();
    _Fm[_qp].setToIdentity();
  }
  unsigned int qp_max = _qpnum;

  auto formatQP = [qp_max](unsigned int qp)
  {
    if (qp_max < 10)
      return std::to_string(qp); // Single digit
    else
      return (qp < 10) ? "0" + std::to_string(qp) : std::to_string(qp); // Two digits
  };

  // If we are using an an externally provided F to recover with (instead of using a solution user
  // object)
  if (isParamValid("F_ext_rec"))
  {
    for (_qp = 0; _qp < n_points; ++_qp)
    {
      _F_store_Fbar[_qp].setToIdentity();
      // int i = 0;
      //
      // for (int i_ind = 0; i_ind < 3; i_ind++)
      //   for (int j_ind = 0; j_ind < 3; j_ind++)
      //   {
      //     std::cout << MetaPhysicL::raw_value((*_F_recover[i])[_qp]) << std::endl;
      //     _F_store_noFbar[_qp](i_ind, j_ind) = (*_F_recover[i])[_qp];
      //     i++;
      //   }
      _F[_qp].setToIdentity();
    }
  }
  else
  {
    if (_recover == true && _volumetric_locking_correction == true)
    {
      ADReal ave_F_det_init = 0;

      std::vector<std::string> indices = {"x", "y", "z"};
      // Get average
      for (_qp = 0; _qp < n_points; ++_qp)
      {
        unsigned int qp_sel = QpMapping::getQP(_qp + 1, _lookup);

        // Populate tensor from solution object
        for (int i_ind = 0; i_ind < 3; i_ind++)
          for (int j_ind = 0; j_ind < 3; j_ind++)
          {
            _F_store_noFbar[_qp](i_ind, j_ind) = _solution_object_ptr->pointValue(
                _t,
                _current_elem->true_centroid(),
                "Fnobar_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                nullptr);
          }
        _F_store_Fbar[_qp] = _F_store_noFbar[_qp];

        ave_F_det_init += _F_store_noFbar[_qp].det() * _JxW[_qp] * _coord[_qp];
      }

      ave_F_det_init /= _current_elem_volume;

      for (_qp = 0; _qp < n_points; ++_qp)
      // Store value
      {
        _F_store_Fbar[_qp] *= cbrt(ave_F_det_init / _F_store_noFbar[_qp].det());
        // For getting the old value of F
        _F[_qp] = _F_store_noFbar[_qp];
      }
    }

    // Recovering without fbar method
    if (_recover == true && _volumetric_locking_correction == false)
    {
      std::vector<std::string> indices = {"x", "y", "z"};
      for (_qp = 0; _qp < n_points; ++_qp)
      {
        unsigned int qp_sel = QpMapping::getQP(_qp + 1, _lookup);

        _F_store_noFbar[_qp].setToIdentity();
        // Populate tensor from solution object
        for (int i_ind = 0; i_ind < 3; i_ind++)
          for (int j_ind = 0; j_ind < 3; j_ind++)
          {
            _F_store_noFbar[_qp](i_ind, j_ind) = _solution_object_ptr->pointValue(
                _t,
                _current_elem->true_centroid(),
                "Fnobar_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                nullptr);
          }
        _F_store_Fbar[_qp] = _F_store_noFbar[_qp];

        // For getting the old value of F
        _F[_qp] = _F_store_noFbar[_qp];
      }
    }
  }
}

void
ComputePlaneDeformationGradient::displacementIntegrityCheck()
{
  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_ndisp != _mesh.dimension())
    paramError(
        "displacements",
        "The number of variables supplied in 'displacements' must match the mesh dimension.");

  // Don't use F-bar in 1D
  if (_ndisp == 1 && _volumetric_locking_correction)
    paramError("volumetric_locking_correction", "has to be set to false for 1-D problems.");

  // Check for RZ
  if (getBlockCoordSystem() == Moose::COORD_RZ && _ndisp != 2)
    paramError("displacements",
               "There must be two displacement variables provided, one in r-direction another in "
               "z-direction");
}

void
ComputePlaneDeformationGradient::initQpStatefulProperties()
{
  _F[_qp].setToIdentity();
  _Fm[_qp].setToIdentity();
}
ADReal
ComputePlaneDeformationGradient::computeQpOutOfPlaneGradDisp()
{
  // std::cout << std::exp(_out_of_plane_strain[_qp]) << std::endl;
  
  //  std::cout << _out_of_plane_strain[_qp] << std::endl;
  // return std::exp(_out_of_plane_strain[_qp]);
    return 1+_out_of_plane_strain[_qp];

}

void
ComputePlaneDeformationGradient::computeProperties()
{
  using std::cbrt;
  ADReal ave_F_det = 0;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    ADRankTwoTensor A = ADRankTwoTensor::initializeFromRows(
        (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    // if (_coord_sys == Moose::COORD_RZ)
    //  std::cout << computeQpOutOfPlaneGradDisp() << std::endl;
    A(2, 2) = computeQpOutOfPlaneGradDisp() - 1;
    _F[_qp] = A;
    _F[_qp].addIa(1.0);

    _Fnobar[_qp].setToIdentity();
    // Outputting the Frobenius norm for post processing reasons
    ADRankTwoTensor temp = _F[_qp];
    temp.addIa(-1);
    _Frobenius[_qp] = MetaPhysicL::raw_value(temp).norm();
    // Add in recovered F
    if (_recover == true)
      _Fnobar[_qp] = _F[_qp] * _F_store_noFbar[_qp];
    else
      _Fnobar[_qp] = _F[_qp];

    if (_volumetric_locking_correction)
      ave_F_det += _F[_qp].det() * _JxW[_qp] * _coord[_qp];
  }

  if (_volumetric_locking_correction)
    ave_F_det /= _current_elem_volume;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    if (_volumetric_locking_correction)
      _F[_qp] *= cbrt(ave_F_det / _F[_qp].det());
    // Multiply in old deformation
    if (_recover == true)
      _F[_qp] = _F[_qp] * _F_store_Fbar[_qp];
    // Remove the eigen deformation gradient
    ADRankTwoTensor Fg(ADRankTwoTensor::initIdentity);
    for (auto Fgi : _Fgs)
      Fg *= (*Fgi)[_qp];
    _Fm[_qp] = Fg.inverse() * _F[_qp];
  }
}
