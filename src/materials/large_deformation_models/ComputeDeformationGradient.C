//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputeDeformationGradient.h"
#include "SolutionUserObjectQP.h"

registerADMooseObject("raccoonApp", ComputeDeformationGradient);

InputParameters
ComputeDeformationGradient::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.addClassDescription(
      "This class computes the deformation gradient. Eigen deformation gradients are extracted "
      "from the total deformation gradient. The F-bar approach can optionally be used to correct "
      "volumetric locking. The deformation of an earlier run can be recovered from a "
      "SolutionUserObject.");

  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<bool>(
      "volumetric_locking_correction", false, "Flag to correct volumetric locking");
  params.addParam<std::vector<MaterialPropertyName>>(
      "eigen_deformation_gradient_names", {}, "List of eigen deformation gradients to be applied");
  params.suppressParameter<bool>("use_displaced_mesh");

  // Polar decomposition output
  params.addParam<bool>("output_half_rotation_tensor",
                        false,
                        "Store R^(1/2) (rotation by theta/2) in the rotation_tensor output "
                        "instead of R. This halves the rotation angle seen by the remapping "
                        "algorithm's matrix logarithm, avoiding the singularity near theta=pi.");
  params.addParam<bool>("use_iterative_polar_decomposition",
                        false,
                        "Use Higham's iterative algorithm for polar decomposition instead of the "
                        "default eigendecomposition of C = F^T*F. Higham's method operates "
                        "directly on F (condition number kappa(F) vs kappa(F)^2 for the eigen "
                        "method), avoids inverting U, requires only 3x3 matrix inverses, and "
                        "converges quadratically. More robust for large stretches and rotations "
                        "near pi.");

  // Recovery
  params.addParam<bool>(
      "recover", false, "Recover the deformation of an earlier run from 'solution'");
  params.addParam<UserObjectName>(
      "solution", "The SolutionUserObjectQP to extract data from; required when recover = true");
  MooseEnum recover_mode("deformation_gradient polar_decomposition", "deformation_gradient");
  params.addParam<MooseEnum>("recover_mode",
                             recover_mode,
                             "Recovery mode: 'deformation_gradient' reads F directly, "
                             "'polar_decomposition' reads R and U and reconstructs F = R*U");
  params.addParam<bool>("input_half_rotation_tensor",
                        false,
                        "Indicates that the rotation tensor in the recovery file contains R^(1/2) "
                        "(i.e., the run that produced the file had output_half_rotation_tensor = "
                        "true). When set, R is reconstructed as R_half * R_half before F = R*U. "
                        "Set to false (default) when the recovery file stores the full R.");
  params.addParam<bool>(
      "recover_apply_fbar_to_U",
      false,
      "How the F-bar corrected deformation gradient is recovered when both 'recover' and "
      "'volumetric_locking_correction' are true. By default (false) the corrected stretch "
      "'stretch_tensor_fbar' is read from the recovery file. When true, the raw stretch "
      "'stretch_tensor' is read instead and the correction is rebuilt on this mesh: F_raw = R*U, "
      "J_avg = <det F_raw>, F_bar = F_raw * cbrt(J_avg / det F_raw). Use this when the recovery "
      "file does not hold stretch_tensor_fbar.");
  params.addParam<bool>(
      "recover_apply_fbar_to_total",
      false,
      "Alternative to 'recover_apply_fbar_to_U' when both 'recover' and "
      "'volumetric_locking_correction' are true: read the raw stretch 'stretch_tensor' and apply "
      "the F-bar correction once to the total raw deformation gradient F_inc * F_rec at every "
      "step, instead of correcting F_inc and F_rec separately and multiplying them. The average "
      "of det F is taken over the original undeformed volume, J_avg = (sum J_inc * JxW) / "
      "(sum JxW / J_rec), which is what a run without a restart computes. Recommended for "
      "explicit dynamics, where no Newton iteration absorbs the error of composing two F-bar "
      "corrections.");
  params.addCoupledVar("F_ext_rec", "External F to recover with");
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
    _recover(getParam<bool>("recover")),
    _solution_object_ptr(nullptr),
    _F_recover(adCoupledValues("F_ext_rec")),
    _Frobenius(declareProperty<Real>(prependBaseName("Frobenius_norm"))),
    _Jacobian(declareADProperty<Real>(prependBaseName("Jacobian"))),
    _rotation_tensor(declareADProperty<RankTwoTensor>(prependBaseName("rotation_tensor"))),
    _rotation_tensor_old(getMaterialPropertyOld<RankTwoTensor>(prependBaseName("rotation_tensor"))),
    _stretch_tensor(declareADProperty<RankTwoTensor>(prependBaseName("stretch_tensor"))),
    _stretch_tensor_fbar(
        declareADProperty<RankTwoTensor>(prependBaseName("stretch_tensor_fbar"))),
    _recover_from_polar(getParam<MooseEnum>("recover_mode") == "polar_decomposition"),
    _output_half_rotation(getParam<bool>("output_half_rotation_tensor")),
    _input_half_rotation(getParam<bool>("input_half_rotation_tensor")),
    _use_iterative_polar(getParam<bool>("use_iterative_polar_decomposition")),
    _recover_apply_fbar_to_U(getParam<bool>("recover_apply_fbar_to_U")),
    _recover_apply_fbar_to_total(getParam<bool>("recover_apply_fbar_to_total")),
    _recover_with_fbar(getParam<bool>("volumetric_locking_correction"))
{
  if (_recover_apply_fbar_to_U && _recover_apply_fbar_to_total)
    paramError("recover_apply_fbar_to_total",
               "Cannot be used together with 'recover_apply_fbar_to_U': they are alternative ways "
               "of recovering the F-bar corrected deformation gradient.");

  for (unsigned int i = 0; i < _Fgs.size(); ++i)
    _Fgs[i] = &Material::getADMaterialProperty<RankTwoTensor>(_Fg_names[i]);

  if (MaterialBase::getParam<bool>("use_displaced_mesh"))
    MaterialBase::paramError("use_displaced_mesh",
                             "The strain calculator needs to run on the undisplaced mesh.");
  if (_recover && !isParamValid("solution"))
    paramError("solution", "Required when recover = true.");
}

void
ComputeDeformationGradient::initialSetup()
{
  displacementIntegrityCheck();

  if (_recover)
    _solution_object_ptr = &getUserObject<SolutionUserObjectQP>("solution");

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

bool
ComputeDeformationGradient::needRawStretch() const
{
  return !_recover_with_fbar || _recover_apply_fbar_to_U || _recover_apply_fbar_to_total;
}

bool
ComputeDeformationGradient::needFbarStretch() const
{
  return !needRawStretch();
}

void
ComputeDeformationGradient::readRecoveredState(std::vector<RankTwoTensor> & R,
                                               std::vector<RankTwoTensor> & R_file,
                                               std::vector<RankTwoTensor> & U,
                                               std::vector<RankTwoTensor> & U_fbar) const
{
  // Only the stretch that is needed is read, so the recovery file has to hold only that one.
  // (It should not hold both: the Exodus reader confuses stretch_tensor_* with
  // stretch_tensor_fbar_* because one name is a prefix of the other.)
  const bool raw = needRawStretch();
  const bool fbar = needFbarStretch();

  const auto read = [&](const std::string & name, const unsigned int qp)
  { return _solution_object_ptr->qpTensor(_t, *_current_elem, name, qp); };

  const unsigned int n_points = R.size();
  for (const auto qp : make_range(n_points))
  {

    if (_recover_from_polar)
    {
      R_file[qp] = read("rotation_tensor", qp);
      R[qp] = _input_half_rotation ? R_file[qp] * R_file[qp] : R_file[qp];
      // The stretch that is not read stands in for the other, so F_raw = R * U and
      // F_fbar = R * U_fbar are both defined
      if (raw)
        U[qp] = read("stretch_tensor", qp);
      U_fbar[qp] = fbar ? read("stretch_tensor_fbar", qp) : U[qp];
      if (!raw)
        U[qp] = U_fbar[qp];
    }
    else
    {
      // Read F and decompose it, so the rest of the recovery only deals with R and U
      RankTwoTensor F_raw, F_fbar;
      if (raw)
        F_raw = read("Fnobar", qp);
      if (fbar)
        F_fbar = read("F", qp);
      const RankTwoTensor & F = raw ? F_raw : F_fbar;
      F.getRUDecompositionRotation(R[qp]);
      R_file[qp] = R[qp];
      U[qp] = R[qp].transpose() * F;
      // F-bar scales F by a scalar, so F_fbar has the same rotation as F_raw
      U_fbar[qp] = fbar ? R[qp].transpose() * F_fbar : U[qp];
    }
  }
}

void
ComputeDeformationGradient::initStatefulProperties(unsigned int n_points)
{
  using std::cbrt;

  for (_qp = 0; _qp < n_points; ++_qp)
  {
    _F[_qp].setToIdentity();
    _Fm[_qp].setToIdentity();
    _rotation_tensor[_qp].setToIdentity();
  }

  // An externally provided F is picked up in computeProperties()
  if (isParamValid("F_ext_rec"))
  {
    for (_qp = 0; _qp < n_points; ++_qp)
      _F_store_Fbar[_qp].setToIdentity();
    return;
  }

  if (!_recover)
    return;

  std::vector<RankTwoTensor> R(n_points), R_file(n_points), U(n_points), U_fbar(n_points);
  readRecoveredState(R, R_file, U, U_fbar);

  // When the F-bar correction is rebuilt from the raw stretch, the element average of det(F_raw)
  // is needed first
  const bool rebuild_fbar =
      _recover_with_fbar && (_recover_apply_fbar_to_U || _recover_apply_fbar_to_total);
  Real J_avg = 0;
  if (rebuild_fbar && _recover_apply_fbar_to_U)
  {
    // Average over this (deformed) mesh
    for (const auto qp : make_range(n_points))
      J_avg += (R[qp] * U[qp]).det() * _JxW[qp] * _coord[qp];
    J_avg /= _current_elem_volume;
  }
  else if (rebuild_fbar)
  {
    // Average over the original undeformed volume, using dV_orig = dV / det(F_raw). This is the
    // average the earlier run computed, so its F-bar is reproduced at each quadrature point.
    Real V = 0, V_orig = 0;
    for (const auto qp : make_range(n_points))
    {
      V += _JxW[qp] * _coord[qp];
      V_orig += _JxW[qp] * _coord[qp] / (R[qp] * U[qp]).det();
    }
    J_avg = V / V_orig;
  }

  for (_qp = 0; _qp < n_points; ++_qp)
  {
    const RankTwoTensor F_raw = R[_qp] * U[_qp];
    RankTwoTensor F_fbar;
    if (!_recover_with_fbar)
      F_fbar = F_raw;
    else if (rebuild_fbar)
      F_fbar = F_raw * cbrt(J_avg / F_raw.det());
    else
      F_fbar = R[_qp] * U_fbar[_qp];

    _F_store_noFbar[_qp] = F_raw;
    _F_store_Fbar[_qp] = F_fbar;

    if (_output_half_rotation)
      _rotation_tensor[_qp] =
          _input_half_rotation
              ? R_file[_qp]
              : computeHalfRotation(R[_qp], RankTwoTensor(RankTwoTensor::initIdentity));
    else
      _rotation_tensor[_qp] = R[_qp];

    _F[_qp] = _recover_with_fbar ? _F_store_Fbar[_qp] : _F_store_noFbar[_qp];

    // Seed Fm with the recovered deformation, so materials that initialize their own stateful
    // properties from it (e.g. the stress) start from the recovered state instead of identity
    ADRankTwoTensor Fg(ADRankTwoTensor::initIdentity);
    for (auto Fgi : _Fgs)
      Fg *= (*Fgi)[_qp];
    _Fm[_qp] = Fg.inverse() * _F[_qp];
  }
}

void
ComputeDeformationGradient::polarDecompositionIterative(const RankTwoTensor & F,
                                                        RankTwoTensor & R,
                                                        RankTwoTensor & U) const
{
  // Higham's Newton iteration for the orthogonal polar factor:
  //   X_{k+1} = (X_k + X_k^{-T}) / 2,  X_0 = F
  // Converges quadratically to R.  Works directly on F (not C = F^T*F), so the
  // effective condition number is kappa(F) rather than kappa(F)^2.  Requires only
  // 3x3 matrix inverses — no LAPACK SVD needed.
  const unsigned int max_iter = 50;
  const Real tol = 1e-12;

  RankTwoTensor X = F;
  for (unsigned int iter = 0; iter < max_iter; ++iter)
  {
    const RankTwoTensor X_new = 0.5 * (X + X.inverse().transpose());
    const Real delta = (X_new - X).norm();
    X = X_new;
    if (delta < tol * X.norm())
      break;
    if (iter == max_iter - 1)
      mooseWarning("polarDecompositionIterative: failed to converge in ",
                   max_iter,
                   " iterations (residual = ",
                   delta,
                   "). Result may be inaccurate.");
  }

  R = X;
  U = R.transpose() * F;
}

RankTwoTensor
ComputeDeformationGradient::computeHalfRotation(const RankTwoTensor & R,
                                                const RankTwoTensor & R_half_old)
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
    K(0, 1) = -n[2];
    K(0, 2) = n[1];
    K(1, 0) = n[2];
    K(1, 2) = -n[0];
    K(2, 0) = -n[1];
    K(2, 1) = n[0];
  }

  // Axis-continuity tracking: when the physical rotation crosses pi, the polar
  // decomposition wraps back and flips the rotation axis sign.  This causes a
  // sudden jump in R^(1/2) components.  We detect the flip by comparing the new
  // axis (encoded in K's skew entries) against the direction implied by the skew
  // part of the previous-timestep R^(1/2).  If the dot product is negative the
  // axis flipped; negate K to restore continuity.
  //
  // The skew part of R_half_old equals sin(theta_old/2)*K_old, so it carries
  // the sign of the old axis without requiring normalisation.  At the first step
  // R_half_old == I and the old skew part is zero — the correction is skipped.
  const Real old_ax = (R_half_old(2, 1) - R_half_old(1, 2)) / 2.0; // ~ sin * n_x
  const Real old_ay = (R_half_old(0, 2) - R_half_old(2, 0)) / 2.0; // ~ sin * n_y
  const Real old_az = (R_half_old(1, 0) - R_half_old(0, 1)) / 2.0; // ~ sin * n_z
  const Real old_mag = std::sqrt(old_ax * old_ax + old_ay * old_ay + old_az * old_az);

  if (old_mag > 1e-10)
  {
    // K(2,1) = n_x, K(0,2) = n_y, K(1,0) = n_z  (from skew-symmetric convention)
    const Real dot = K(2, 1) * old_ax + K(0, 2) * old_ay + K(1, 0) * old_az;
    if (dot < 0.0)
      K *= -1.0;
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
ComputeDeformationGradient::applyOutOfPlaneGradDisp(ADRankTwoTensor & A)
{
  // Axisymmetric: F_zz is the hoop stretch u_r / r. Cartesian bulk problems leave A(2,2) as the
  // (zero) out-of-plane displacement gradient. Plane-stress models override this hook.
  if (_coord_sys == Moose::COORD_RZ)
    A(2, 2) = computeQpOutOfPlaneGradDisp();
}

void
ComputeDeformationGradient::checkPolarDecomposition(const RankTwoTensor & R,
                                                    const RankTwoTensor & U,
                                                    const RankTwoTensor & F) const
{
  const Real tol = 1e-6;
  const auto warn = [&](const std::string & quantity, const Real value, const std::string & what)
  {
    mooseWarning(name(),
                 ": polar decomposition ",
                 quantity,
                 " = ",
                 value,
                 " at QP ",
                 _qp,
                 " element ",
                 _current_elem->id(),
                 ". ",
                 what);
  };

  const Real det_R = R.det();
  if (std::abs(det_R - 1.0) > tol)
    warn("det(R)", det_R, "R is not a proper rotation.");

  const Real ortho_norm = (R * R.transpose() - RankTwoTensor(RankTwoTensor::initIdentity)).norm();
  if (ortho_norm > tol)
    warn("||R*R^T - I||", ortho_norm, "R is not orthogonal.");

  const Real recon_norm = (R * U - F).norm();
  if (recon_norm > tol)
    warn("||R*U - F||", recon_norm, "Reconstruction error is large.");
}

void
ComputeDeformationGradient::computeProperties()
{
  using std::cbrt;

  // Whether F-bar is applied once to the total raw deformation gradient (see the class description)
  const bool fbar_to_total =
      _recover && _recover_apply_fbar_to_total && _volumetric_locking_correction;

  // Volume integrals for the element average of det(F)
  ADReal ave_F_det = 0;
  ADReal original_volume = 0;

  if (isParamValid("F_ext_rec") && _t_step >= 1)
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      for (const auto i : make_range(3))
        for (const auto j : make_range(3))
          _F_store_Fbar[_qp](i, j) = MetaPhysicL::raw_value((*_F_recover[3 * i + j])[_qp]);

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    // The (incremental, when recovering) deformation gradient from the displacements
    ADRankTwoTensor A = ADRankTwoTensor::initializeFromRows(
        (*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);
    applyOutOfPlaneGradDisp(A);
    _F[_qp] = A;
    _F[_qp].addIa(1.0);

    // For postprocessing
    _Frobenius[_qp] = MetaPhysicL::raw_value(A).norm();
    _Jacobian[_qp] = _F[_qp].det();

    // The total deformation gradient without the F-bar correction
    _Fnobar[_qp] = _recover ? _F[_qp] * _F_store_noFbar[_qp] : _F[_qp];

    if (_volumetric_locking_correction)
    {
      ave_F_det += _F[_qp].det() * _JxW[_qp] * _coord[_qp];
      if (fbar_to_total)
        original_volume += _JxW[_qp] * _coord[_qp] / _F_store_noFbar[_qp].det();
    }
  }

  // This mesh is the deformed configuration of the earlier run when recovering, so dividing by
  // the original volume (dV_orig = dV / det(F_rec)) gives the average of det(F_total) over the
  // original volume, as in a run without a restart.
  if (_volumetric_locking_correction)
    ave_F_det /= fbar_to_total ? original_volume : ADReal(_current_elem_volume);

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    if (fbar_to_total)
      _F[_qp] = _Fnobar[_qp] * cbrt(ave_F_det / _Fnobar[_qp].det());
    else
    {
      if (_volumetric_locking_correction)
        _F[_qp] *= cbrt(ave_F_det / _F[_qp].det());

      // Multiply in the recovered deformation
      if (_recover)
        _F[_qp] = _F[_qp] * _F_store_Fbar[_qp];
    }

    // Remove the eigen deformation gradient
    ADRankTwoTensor Fg(ADRankTwoTensor::initIdentity);
    for (auto Fgi : _Fgs)
      Fg *= (*Fgi)[_qp];
    _Fm[_qp] = Fg.inverse() * _F[_qp];

    // Polar decomposition F = R * U of the deformation gradient without the F-bar correction
    const RankTwoTensor Fnobar = MetaPhysicL::raw_value(_Fnobar[_qp]);
    RankTwoTensor R, U;
    if (_use_iterative_polar)
      polarDecompositionIterative(Fnobar, R, U);
    else
    {
      Fnobar.getRUDecompositionRotation(R);
      U = R.transpose() * Fnobar;
    }
    checkPolarDecomposition(R, U, Fnobar);

    _rotation_tensor[_qp] =
        _output_half_rotation ? computeHalfRotation(R, _rotation_tensor_old[_qp]) : R;
    _stretch_tensor[_qp] = U;

    // F-bar scales F by a scalar, so F_fbar = R * U_fbar with the same rotation. A restart reads
    // U_fbar to recover F_fbar without averaging det(F) again on its own mesh.
    _stretch_tensor_fbar[_qp] = ADRankTwoTensor(R).transpose() * _F[_qp];
  }
}
