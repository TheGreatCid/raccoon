//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "ADReal.h"
#include "ComputeDeformationGradient.h"
#include "EigenADReal.h"
#include "Material.h"
#include "Moose.h"
#include "MooseError.h"
#include "MooseTypes.h"
#include "RankTwoTensorForward.h"
#include "SolutionUserObject.h"
#include "metaphysicl/raw_type.h"
#include "Qp_Mapping.h"

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
  MooseEnum recover_mode("deformation_gradient polar_decomposition", "deformation_gradient");
  params.addParam<MooseEnum>("recover_mode",
                             recover_mode,
                             "Recovery mode: 'deformation_gradient' reads F directly, "
                             "'polar_decomposition' reads R and U and reconstructs F = R*U");
  params.addParam<bool>("output_half_rotation_tensor",
                        false,
                        "Store R^(1/2) (rotation by theta/2) in the rotation_tensor output "
                        "instead of R. This halves the rotation angle seen by the remapping "
                        "algorithm's matrix logarithm, avoiding the singularity near theta=pi.");
  params.addParam<bool>("input_half_rotation_tensor",
                        false,
                        "Indicates that the rotation tensor in the recovery file contains R^(1/2) "
                        "(i.e., the run that produced the file had output_half_rotation_tensor = "
                        "true). When set, R is reconstructed as R_half * R_half before F = R*U. "
                        "Set to false (default) when the recovery file stores the full R.");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<UserObjectName>("solution", "The SolutionUserObject to extract data from.");
  params.addParam<Real>("num_qps", 8, "Number of QPs");
  params.addCoupledVar("F_ext_rec", "External F to recover with");
  params.addRequiredParam<MooseEnum>(
      "element", MooseEnum(QpMapping::ELEMENT_ENUM_DEFINITION), "The element type");
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
    _weights(declareADProperty<Real>("weights")),
    _recover(getParam<bool>("recover")),
    _solution_object_ptr(NULL),
    _F_recover(adCoupledValues("F_ext_rec")),
    _element(getParam<MooseEnum>("element").getEnum<QpMapping::Element>()),
    _Frobenius(declareProperty<Real>(prependBaseName("Frobenius_norm"))),
    _Jacobian(declareProperty<Real>(prependBaseName("Jacobian"))),
    _rotation_tensor(declareADProperty<RankTwoTensor>(prependBaseName("rotation_tensor"))),
    _stretch_tensor(declareADProperty<RankTwoTensor>(prependBaseName("stretch_tensor"))),
    _recover_from_polar(getParam<MooseEnum>("recover_mode") == "polar_decomposition"),
    _output_half_rotation(getParam<bool>("output_half_rotation_tensor")),
    _input_half_rotation(getParam<bool>("input_half_rotation_tensor"))
{
  for (unsigned int i = 0; i < _Fgs.size(); ++i)
    _Fgs[i] = &Material::getADMaterialProperty<RankTwoTensor>(_Fg_names[i]);

  if (MaterialBase::getParam<bool>("use_displaced_mesh"))
    MaterialBase::paramError("use_displaced_mesh",
                             "The strain calculator needs to run on the undisplaced mesh.");
  if (_recover)
    _lookup = QpMapping::getLookup(_element, _qpnum, /*reversed=*/true);
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

        if (_recover_from_polar)
        {
          // Recover from rotation and stretch tensors
          RankTwoTensor R, U;
          for (int i_ind = 0; i_ind < 3; i_ind++)
            for (int j_ind = 0; j_ind < 3; j_ind++)
            {
              R(i_ind, j_ind) = _solution_object_ptr->pointValue(
                  _t,
                  _current_elem->true_centroid(),
                  "rotation_tensor_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                  nullptr);
              U(i_ind, j_ind) = _solution_object_ptr->pointValue(
                  _t,
                  _current_elem->true_centroid(),
                  "stretch_tensor_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                  nullptr);
            }
          // Reconstruct F = R * U and convert to AD type.
          // If the recovery file stored R^(1/2), square it first to recover full R.
          if (_input_half_rotation)
            R = R * R;
          RankTwoTensor F_reconstructed = R * U;
          for (int i_ind = 0; i_ind < 3; i_ind++)
            for (int j_ind = 0; j_ind < 3; j_ind++)
              _F_store_noFbar[_qp](i_ind, j_ind) = F_reconstructed(i_ind, j_ind);
        }
        else
        {
          // Populate tensor from solution object (traditional recovery from F)
          for (int i_ind = 0; i_ind < 3; i_ind++)
            for (int j_ind = 0; j_ind < 3; j_ind++)
            {
              _F_store_noFbar[_qp](i_ind, j_ind) = _solution_object_ptr->pointValue(
                  _t,
                  _current_elem->true_centroid(),
                  "Fnobar_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                  nullptr);
            }
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

        if (_recover_from_polar)
        {
          // Recover from rotation and stretch tensors
          RankTwoTensor R, U;
          for (int i_ind = 0; i_ind < 3; i_ind++)
            for (int j_ind = 0; j_ind < 3; j_ind++)
            {
              R(i_ind, j_ind) = _solution_object_ptr->pointValue(
                  _t,
                  _current_elem->true_centroid(),
                  "rotation_tensor_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                  nullptr);
              U(i_ind, j_ind) = _solution_object_ptr->pointValue(
                  _t,
                  _current_elem->true_centroid(),
                  "stretch_tensor_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                  nullptr);
            }
          // Reconstruct F = R * U and convert to AD type.
          // If the recovery file stored R^(1/2), square it first to recover full R.
          if (_input_half_rotation)
            R = R * R;
          RankTwoTensor F_reconstructed = R * U;
          for (int i_ind = 0; i_ind < 3; i_ind++)
            for (int j_ind = 0; j_ind < 3; j_ind++)
              _F_store_noFbar[_qp](i_ind, j_ind) = F_reconstructed(i_ind, j_ind);
        }
        else
        {
          // Populate tensor from solution object (traditional recovery from F)
          for (int i_ind = 0; i_ind < 3; i_ind++)
            for (int j_ind = 0; j_ind < 3; j_ind++)
            {
              _F_store_noFbar[_qp](i_ind, j_ind) = _solution_object_ptr->pointValue(
                  _t,
                  _current_elem->true_centroid(),
                  "Fnobar_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                  nullptr);
            }
        }

        _F_store_Fbar[_qp] = _F_store_noFbar[_qp];

        // For getting the old value of F
        _F[_qp] = _F_store_noFbar[_qp];
      }
    }
  }
}

RankTwoTensor
ComputeDeformationGradient::computeHalfRotation(const RankTwoTensor & R)
{
  // Extract rotation angle from the trace: cos(theta) = (tr(R) - 1) / 2
  const Real cos_theta = std::max(-1.0, std::min(1.0, (R.tr() - 1.0) / 2.0));
  const Real theta = std::acos(cos_theta);

  const RankTwoTensor I(RankTwoTensor::initIdentity);

  if (theta < 1e-10)
    return I;

  RankTwoTensor K; // skew-symmetric cross-product matrix of the rotation axis

  if (std::abs(theta - M_PI) > 1e-4)
  {
    // General case: extract axis from the skew-symmetric part of R.
    // (R - R^T) / 2 = sin(theta) * K
    const Real s = std::sin(theta);
    for (const auto i : make_range(3))
      for (const auto j : make_range(3))
        K(i, j) = (R(i, j) - R(j, i)) / (2.0 * s);
  }
  else
  {
    // Near theta = pi: skew part vanishes (sin(pi) = 0), extract axis from
    // the symmetric part instead.  At theta = pi: R = -I + 2*n*n^T, so
    // (R + I)/2 = n*n^T.  Find the dominant diagonal entry to get n_max,
    // then recover the remaining components from the off-diagonal.
    int idx = 0;
    Real max_diag = (R(0, 0) + 1.0) / 2.0;
    for (const auto k : make_range(1, 3))
    {
      const Real d = (R(k, k) + 1.0) / 2.0;
      if (d > max_diag)
      {
        max_diag = d;
        idx = k;
      }
    }
    Real n[3] = {0.0, 0.0, 0.0};
    n[idx] = std::sqrt(std::max(0.0, max_diag));
    for (const auto k : make_range(3))
      if (k != idx)
        n[k] = R(idx, k) / (2.0 * n[idx]);
    // Normalise for robustness
    const Real len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (const auto k : make_range(3))
      n[k] /= len;
    K(0, 1) = -n[2]; K(0, 2) =  n[1];
    K(1, 0) =  n[2]; K(1, 2) = -n[0];
    K(2, 0) = -n[1]; K(2, 1) =  n[0];
  }

  // Rodrigues formula for R^(1/2): same axis, half angle.
  // R = I + sin(theta)*K + (1-cos(theta))*K^2
  // R^(1/2) = I + sin(theta/2)*K + (1-cos(theta/2))*K^2
  const Real half = theta / 2.0;
  return I + std::sin(half) * K + (1.0 - std::cos(half)) * (K * K);
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
  using std::cbrt;

  ADReal ave_F_det = 0;

  if (isParamValid("F_ext_rec"))
  {
    if (_t_step >= 1)
    {
      for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      {
        int i = 0;
        for (int i_ind = 0; i_ind < 3; i_ind++)
          for (int j_ind = 0; j_ind < 3; j_ind++)
          {
            _F_store_Fbar[_qp](i_ind, j_ind) = MetaPhysicL::raw_value((*_F_recover[i])[_qp]);
            i++;
          }
      }
    }
  }

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    ADRankTwoTensor A = ADRankTwoTensor::initializeFromRows(
        (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    if (_coord_sys == Moose::COORD_RZ)
      A(2, 2) = computeQpOutOfPlaneGradDisp();
    _F[_qp] = A;
    _F[_qp].addIa(1.0);

    _Fnobar[_qp].setToIdentity();

    // Outputting the Frobenius norm for post processing reasons
    ADRankTwoTensor temp = _F[_qp];
    temp.addIa(-1);
    _Frobenius[_qp] = MetaPhysicL::raw_value(temp).norm();

    // Outputting the Jacobian (determinant of F) for post processing reasons
    _Jacobian[_qp] = MetaPhysicL::raw_value(_F[_qp].det());

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

    // Compute polar decomposition of non-volume corrected F
    // Convert ADRankTwoTensor to RankTwoTensor for polar decomposition
    RankTwoTensor Fnobar_real = MetaPhysicL::raw_value(_Fnobar[_qp]);
    RankTwoTensor R, U;

    // Compute rotation tensor R from F = R*U decomposition
    Fnobar_real.getRUDecompositionRotation(R);

    // Compute right stretch tensor U = R^T * F
    U = R.transpose() * Fnobar_real;

    // Verify polar decomposition quality
    {
      const Real det_R = R.det();
      if (std::abs(det_R - 1.0) > 1e-6)
        mooseWarning(name(),
                     ": polar decomposition det(R) = ",
                     det_R,
                     " (expected 1) at QP ",
                     _qp,
                     " element ",
                     _current_elem->id(),
                     ". R is not a proper rotation.");

      const RankTwoTensor ortho_err =
          R * R.transpose() - RankTwoTensor(RankTwoTensor::initIdentity);
      const Real ortho_norm = ortho_err.norm();
      if (ortho_norm > 1e-6)
        mooseWarning(name(),
                     ": polar decomposition ||R*R^T - I|| = ",
                     ortho_norm,
                     " at QP ",
                     _qp,
                     " element ",
                     _current_elem->id(),
                     ". R is not orthogonal.");

      const Real recon_norm = (R * U - Fnobar_real).norm();
      if (recon_norm > 1e-6)
        mooseWarning(name(),
                     ": polar decomposition ||R*U - F|| = ",
                     recon_norm,
                     " at QP ",
                     _qp,
                     " element ",
                     _current_elem->id(),
                     ". Reconstruction error is large.");
    }

    _rotation_tensor[_qp] = _output_half_rotation ? computeHalfRotation(R) : R;
    _stretch_tensor[_qp] = U;
  }
}
