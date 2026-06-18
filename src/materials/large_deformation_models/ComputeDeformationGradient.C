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
  params.addParam<bool>("use_iterative_polar_decomposition",
                        false,
                        "Use Higham's iterative algorithm for polar decomposition instead of the "
                        "default eigendecomposition of C = F^T*F. Higham's method operates "
                        "directly on F (condition number kappa(F) vs kappa(F)^2 for the eigen "
                        "method), avoids inverting U, requires only 3x3 matrix inverses, and "
                        "converges quadratically. More robust for large stretches and rotations "
                        "near pi.");
  params.addParam<bool>("input_half_rotation_tensor",
                        false,
                        "Indicates that the rotation tensor in the recovery file contains R^(1/2) "
                        "(i.e., the run that produced the file had output_half_rotation_tensor = "
                        "true). When set, R is reconstructed as R_half * R_half before F = R*U. "
                        "Set to false (default) when the recovery file stores the full R.");
  params.addParam<bool>(
      "recover_apply_fbar_to_U",
      false,
      "When recovering with volumetric_locking_correction = true, rebuild F_bar on the new mesh "
      "by applying the original F-bar volumetric-averaging operator to the recovered raw U "
      "(F_raw = R*U; J_avg = <det F_raw>; F_bar = F_raw * cbrt(J_avg/det F_raw)) -- approach B "
      "from the cfb539fe1-era recovery method.  When false (default), F_bar is read directly "
      "from the recovery file's `stretch_tensor_fbar` (approach A: recover U_bar already "
      "averaged on the reference side).  Use approach B when the recovery file does not carry "
      "stretch_tensor_fbar and the F-bar correction must be re-derived on the new mesh.  "
      "Active only when `recover = true` and `volumetric_locking_correction = true`.");
  params.suppressParameter<bool>("use_displaced_mesh");
  params.addParam<UserObjectName>("solution", "The SolutionUserObject to extract data from.");
  params.addParam<Real>("num_qps", 8, "Number of QPs");
  params.addCoupledVar("F_ext_rec", "External F to recover with");
  params.addParam<MooseEnum>("element",
                             MooseEnum(QpMapping::ELEMENT_ENUM_DEFINITION),
                             "Element type for QP remapping; required when recover = true");
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
    _element(QpMapping::Element::HEX8_3rd),
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
    _recover_apply_fbar_to_U(getParam<bool>("recover_apply_fbar_to_U"))
{
  for (unsigned int i = 0; i < _Fgs.size(); ++i)
    _Fgs[i] = &Material::getADMaterialProperty<RankTwoTensor>(_Fg_names[i]);

  if (MaterialBase::getParam<bool>("use_displaced_mesh"))
    MaterialBase::paramError("use_displaced_mesh",
                             "The strain calculator needs to run on the undisplaced mesh.");
  if (_recover)
  {
    if (!isParamSetByUser("element"))
      mooseError("'element' must be specified when recover = true");
    _element = getParam<MooseEnum>("element").getEnum<QpMapping::Element>();
    _lookup = QpMapping::getLookup(_element, _qpnum, /*reversed=*/true);
  }
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
    _rotation_tensor[_qp].setToIdentity();
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
    if (_recover == true)
    {
      // Recovery pipeline:
      //   (1) read R/U/U_fbar (or polar-decompose F/F_fbar) at every QP of the element,
      //   (2) reconstruct F_raw = R * U and either pull F_fbar = R * U_fbar from the file
      //       (approach A) or rebuild it via the F-bar averaging operator on F_raw
      //       (approach B; recover_apply_fbar_to_U = true), then store.
      const std::vector<std::string> indices = {"x", "y", "z"};
      const bool have_fbar = _volumetric_locking_correction;
      // Approach B (recover_apply_fbar_to_U) rebuilds F_bar from the raw recovered U on
      // the new mesh -- it doesn't need U_fbar from the recovery file, so don't try to
      // read it.  Approach A (recover_apply_fbar_to_U = false, F-bar active) consumes
      // U_fbar directly and doesn't need the raw U either -- F_noFbar gets aliased to
      // F_Fbar at INITIAL in that case.  Reading only what's needed lets the same
      // restart input work against reference dumps that contain only one of the two
      // stretch-tensor variants (libmesh's exodus reader silently shadows
      // `stretch_tensor_*` whenever `stretch_tensor_fbar_*` is also present in the same
      // file due to a name-prefix collision, so the reference can dump only one).
      const bool need_U_fbar_from_file = have_fbar && !_recover_apply_fbar_to_U;
      const bool need_U_raw_from_file = !have_fbar || _recover_apply_fbar_to_U;

      std::vector<RankTwoTensor> R_qp(n_points), U_qp(n_points), Ufb_qp(n_points);
      std::vector<RankTwoTensor> R_file_qp(n_points); // raw R as read; may be R^(1/2)

      for (unsigned int qp = 0; qp < n_points; ++qp)
      {
        const unsigned int qp_sel = QpMapping::getQP(qp + 1, _lookup);

        if (_recover_from_polar)
        {
          RankTwoTensor R, U, Ufb;
          for (int i_ind = 0; i_ind < 3; ++i_ind)
            for (int j_ind = 0; j_ind < 3; ++j_ind)
            {
              R(i_ind, j_ind) = _solution_object_ptr->pointValue(
                  _t,
                  _current_elem->true_centroid(),
                  "rotation_tensor_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                  nullptr);
              if (need_U_raw_from_file)
                U(i_ind, j_ind) = _solution_object_ptr->pointValue(
                    _t,
                    _current_elem->true_centroid(),
                    "stretch_tensor_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                    nullptr);
              if (need_U_fbar_from_file)
                Ufb(i_ind, j_ind) = _solution_object_ptr->pointValue(
                    _t,
                    _current_elem->true_centroid(),
                    "stretch_tensor_fbar_" + indices[i_ind] + indices[j_ind] + "_" +
                        formatQP(qp_sel),
                    nullptr);
            }
          R_file_qp[qp] = R;
          if (_input_half_rotation)
            R = R * R;
          R_qp[qp] = R;
          // Approach A skips reading raw U; alias U <- U_fbar so F_raw = R*U downstream
          // collapses to F_fbar (i.e., _F_store_noFbar = _F_store_Fbar at INITIAL).
          U_qp[qp] = need_U_raw_from_file ? U : Ufb;
          // Approach B / no-fbar skip U_fbar; alias Ufb <- U so the F_fbar pathway has
          // a sane placeholder.  Approach B overrides _F_store_Fbar below via averaging.
          Ufb_qp[qp] = need_U_fbar_from_file ? Ufb : U_qp[qp];
        }
        else
        {
          // Traditional recovery from F: read F_raw (and F_fbar) directly, then polar-decompose
          // so the downstream code uses a uniform (R, U) view regardless of recovery mode.
          RankTwoTensor F_raw, F_fbar;
          for (int i_ind = 0; i_ind < 3; ++i_ind)
            for (int j_ind = 0; j_ind < 3; ++j_ind)
            {
              if (need_U_raw_from_file)
                F_raw(i_ind, j_ind) = _solution_object_ptr->pointValue(
                    _t,
                    _current_elem->true_centroid(),
                    "Fnobar_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                    nullptr);
              if (need_U_fbar_from_file)
                F_fbar(i_ind, j_ind) = _solution_object_ptr->pointValue(
                    _t,
                    _current_elem->true_centroid(),
                    "F_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
                    nullptr);
            }
          // For approach A (need_U_raw=false but need_U_fbar=true) we alias F_raw <- F_fbar
          // so the polar decomposition still has a valid input.
          const RankTwoTensor & F_for_RU = need_U_raw_from_file ? F_raw : F_fbar;
          RankTwoTensor R, U;
          F_for_RU.getRUDecompositionRotation(R);
          U = R.transpose() * F_for_RU;
          R_file_qp[qp] = R;
          R_qp[qp] = R;
          U_qp[qp] = U;
          // F_fbar = R * U_fbar with the same R since fbar is a scalar volumetric multiplier.
          Ufb_qp[qp] = need_U_fbar_from_file ? R.transpose() * F_fbar : U;
        }
      }

      // Approach B (recover_apply_fbar_to_U = true, F-bar active): rebuild F_bar on the
      // new mesh by applying the original cfb539fe1-era F-bar averaging operator to the
      // recovered raw U.  Pre-pass over QPs computes the element-volume-weighted average
      // det(F_raw); the per-QP F_fbar in the next loop is then scaled to enforce that
      // det(F_bar) = J_avg per element.
      Real J_avg_init = 0.0;
      const bool fbar_from_raw_U = have_fbar && _recover_apply_fbar_to_U;
      if (fbar_from_raw_U)
      {
        for (unsigned int qp = 0; qp < n_points; ++qp)
        {
          const RankTwoTensor F_raw_qp = R_qp[qp] * U_qp[qp];
          J_avg_init += F_raw_qp.det() * _JxW[qp] * _coord[qp];
        }
        J_avg_init /= _current_elem_volume;
      }

      for (_qp = 0; _qp < n_points; ++_qp)
      {
        const RankTwoTensor F_raw = R_qp[_qp] * U_qp[_qp];
        // Three branches:
        //   no F-bar           -> F_fbar = F_raw
        //   approach A (default) -> F_fbar = R * U_fbar from the recovery file
        //   approach B           -> F_fbar = F_raw * cbrt(J_avg / det F_raw)  (cfb539fe1)
        RankTwoTensor F_fbar;
        if (!have_fbar)
          F_fbar = F_raw;
        else if (fbar_from_raw_U)
          F_fbar = F_raw * cbrt(J_avg_init / F_raw.det());
        else
          F_fbar = R_qp[_qp] * Ufb_qp[_qp];

        for (int i_ind = 0; i_ind < 3; ++i_ind)
          for (int j_ind = 0; j_ind < 3; ++j_ind)
          {
            _F_store_noFbar[_qp](i_ind, j_ind) = F_raw(i_ind, j_ind);
            _F_store_Fbar[_qp](i_ind, j_ind) = F_fbar(i_ind, j_ind);
          }

        const RankTwoTensor I_seed(RankTwoTensor::initIdentity);
        if (_output_half_rotation)
          _rotation_tensor[_qp] =
              _input_half_rotation ? R_file_qp[_qp] : computeHalfRotation(R_qp[_qp], I_seed);
        else
          _rotation_tensor[_qp] = R_qp[_qp];

        _F[_qp] = have_fbar ? _F_store_Fbar[_qp] : _F_store_noFbar[_qp];

        // Seed _Fm at INITIAL so downstream stateful seeders (e.g.
        // ComputeLargeDeformationStress::initQpStatefulProperties) see
        // F_m = Fg^-1 * F_recovered rather than identity.  Without this,
        // _stress at INITIAL gets seeded from the constitutive on F = I,
        // making _stress_old = 0 at t_step=1 and breaking the HHT-alpha
        // residual on the first restart step.
        ADRankTwoTensor Fg(ADRankTwoTensor::initIdentity);
        for (auto Fgi : _Fgs)
          Fg *= (*Fgi)[_qp];
        _Fm[_qp] = Fg.inverse() * _F[_qp];
      }
    }
  }
}


void
ComputeDeformationGradient::polarDecompositionIterative(const RankTwoTensor & F, // NOLINT
                                                        RankTwoTensor & R,
                                                        RankTwoTensor & U)
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
    _Jacobian[_qp] = _F[_qp].det();

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

    // Compute polar decomposition F = R*U
    if (_use_iterative_polar)
      polarDecompositionIterative(Fnobar_real, R, U);
    else
    {
      Fnobar_real.getRUDecompositionRotation(R);
      U = R.transpose() * Fnobar_real;
    }

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

    _rotation_tensor[_qp] =
        _output_half_rotation ? computeHalfRotation(R, _rotation_tensor_old[_qp]) : R;
    _stretch_tensor[_qp] = U;

    // U_fbar = R^T * F_fbar.  R is from the polar decomposition of F_raw = Fnobar; the
    // volumetric scalar that produces F_fbar from F_raw commutes with the rotation, so
    // F_fbar = R * U_fbar with the same R.  U_fbar is SPD (det = <J>_elem) and on the
    // same manifold as U, so a downstream Lie-algebra interpolator handles it the same
    // way (log -> interp -> exp).  Restarts read U_fbar directly to seed _F_store_Fbar
    // without re-averaging <J> on the new mesh's _JxW.
    ADRankTwoTensor R_ad;
    for (int i_ind = 0; i_ind < 3; ++i_ind)
      for (int j_ind = 0; j_ind < 3; ++j_ind)
        R_ad(i_ind, j_ind) = R(i_ind, j_ind);
    _stretch_tensor_fbar[_qp] = R_ad.transpose() * _F[_qp];
  }
}
