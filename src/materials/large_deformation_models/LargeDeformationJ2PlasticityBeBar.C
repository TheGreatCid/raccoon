//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "EigenADReal.h"
#include "LargeDeformationJ2PlasticityBeBar.h"
#include "CNHIsotropicElasticity.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", LargeDeformationJ2PlasticityBeBar);

InputParameters
LargeDeformationJ2PlasticityBeBar::validParams()
{
  InputParameters params = LargeDeformationJ2PlasticityBase::validParams();
  params.addClassDescription("Large deformation $J_2$ plasticity using the bebar (modified "
                             "left Cauchy-Green) update. Elastic parameters are sourced from "
                             "the associated CNHIsotropicElasticity model.");
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density_corr",
      "psie_corr",
      "Name of the strain energy density computed by this material model");
  params.addParam<bool>("apply_strain_energy_split",
                        true,
                        "Whether to apply volumetric/deviatoric split to the active strain "
                        "energy (psie_active). If true (default), only the positive part of the "
                        "energy is assigned to psie_active. If false, the full unsplit energy U + "
                        "W is used for psie_active.");
  params.addParam<Real>("triax_gaussian_peak",
                        2.0,
                        "Peak value of the Gaussian triaxiality function f(η) at η=0. "
                        "f(η) = 1 + (peak-1)*exp(-η²/(2σ²)), so the asymptotic value at "
                        "|η|→∞ is 1 and f(0)=peak. Default 2 recovers the original behavior.");
  params.addParam<MaterialPropertyName>(
      "hardening_plastic_energy_density_active",
      "psip_active",
      "Name of the active plastic energy density declared by the hardening model "
      "(used for psip_triax). Must match the hardening model's declared property name.");
  params.addParam<Real>("psip_triax_threshold",
                        200.0,
                        "Activation threshold for psip_triax. psip_triax_raw must exceed "
                        "this value before psip_triax becomes non-zero.");
  return params;
}

LargeDeformationJ2PlasticityBeBar::LargeDeformationJ2PlasticityBeBar(
    const InputParameters & parameters)
  : LargeDeformationJ2PlasticityBase(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _bebar(declareADProperty<RankTwoTensor>(prependBaseName("be_bar"))),
    _bebar_old(getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("be_bar"))),
    _bebar_det(declareADProperty<Real>(prependBaseName("be_bar_det"))),
    _F_old(getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _F(getADMaterialPropertyByName<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _d_name(getVar("phase_field", 0)->name()),
    _psie_name(prependBaseName("strain_energy_density_corr", true)),
    _psie_corr(declareADProperty<Real>(_psie_name)),
    _psie_active_corr(declareADProperty<Real>(_psie_name + "_active")),
    _dpsie_dd_corr(declareADProperty<Real>(derivativePropertyName(_psie_name, {_d_name}))),
    _psie_unsplit(declareADProperty<Real>(_psie_name + "_unsplit")),
    _apply_strain_energy_split(getParam<bool>("apply_strain_energy_split")),
    _triaxfunc(declareADProperty<Real>(prependBaseName("triaxfunc"))),
    _triaxiality_kirchhoff(declareADProperty<Real>(prependBaseName("triaxiality_kirchhoff"))),
    _triaxiality_cauchy(declareADProperty<Real>(prependBaseName("triaxiality_cauchy"))),
    _triax_gaussian_peak(getParam<Real>("triax_gaussian_peak")),
    _psip_active_ref(getADMaterialProperty<Real>(
        getParam<MaterialPropertyName>("hardening_plastic_energy_density_active"))),
    _psip_triax_threshold(getParam<Real>("psip_triax_threshold")),
    _psip_triax_raw(declareADProperty<Real>(prependBaseName("psip_triax_raw"))),
    _psip_triax_raw_old(getMaterialPropertyOld<Real>(prependBaseName("psip_triax_raw"))),
    _psip_triax(declareADProperty<Real>(prependBaseName("psip_triax"))),
    _psip_triax_old(getMaterialPropertyOld<Real>(prependBaseName("psip_triax")))
{
  _check_range = true;
}

void
LargeDeformationJ2PlasticityBeBar::setElasticityModel(
    LargeDeformationElasticityModel * elasticity_model)
{
  auto * cnh = dynamic_cast<CNHIsotropicElasticity *>(elasticity_model);
  if (!cnh)
    mooseError(type(), ": requires a CNHIsotropicElasticity model; got ", elasticity_model->type());
  _K = &cnh->getK();
  _G = &cnh->getG();
  _ge = &cnh->getDegradation();
  _dge_dd = &cnh->getDegradationDerivative();
  _psie_cnh = &cnh->getPsie();
  _psie_active_cnh = &cnh->getPsieActive();
  _dpsie_dd_cnh = &cnh->getDpsieDD();
  LargeDeformationJ2PlasticityBase::setElasticityModel(elasticity_model);
}

void
LargeDeformationJ2PlasticityBeBar::initQpStatefulProperties()
{
  _Fp[_qp].setToIdentity();
  _ep[_qp] = 0;
  _bebar[_qp].setToIdentity();
  _psip_triax_raw[_qp] = 0.0;
  _psip_triax[_qp] = 0.0;

  if (!_recover)
    return;

  const std::vector<std::string> indices = {"x", "y", "z"};
  unsigned int qp_sel = QpMapping::getQP(_qp + 1, _lookup);

  _ep[_qp] = _solution_object_ptr->pointValue(
      _t, _current_elem->true_centroid(), "effective_plastic_strain_" + formatQP(qp_sel), nullptr);
  if (_ep[_qp] < 0)
    _ep[_qp] = 0;

  for (int i_ind = 0; i_ind < 3; i_ind++)
    for (int j_ind = 0; j_ind < 3; j_ind++)
      _bebar[_qp](i_ind, j_ind) = _solution_object_ptr->pointValue(
          _t,
          _current_elem->true_centroid(),
          "be_bar_" + indices[i_ind] + indices[j_ind] + "_" + formatQP(qp_sel),
          nullptr);

  // Re-enforce det(be_bar) = 1 on the mapped be_bar.
  const Real det_bebar = MetaPhysicL::raw_value(_bebar[_qp].det());
  if (std::abs(det_bebar - 1.0) > 1e-10)
    _bebar[_qp] /= std::cbrt(det_bebar);
}

void
LargeDeformationJ2PlasticityBeBar::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & /*Fe*/)
{
  using std::cbrt;
  using std::sqrt;

  ADRankTwoTensor I2;
  I2.setToIdentity();

  ADReal delta_ep = 0;

  // Compute incremental fbar
  ADRankTwoTensor f = _F[_qp] * _F_old[_qp].inverse();
  ADRankTwoTensor fbar = f / cbrt(f.det());

  // Compute bebar_trial
  _bebar[_qp] = fbar * _bebar_old[_qp] * fbar.transpose();

  // Compute trial deviatoric stress and its norm
  ADRankTwoTensor s_trial = (*_G)[_qp] * (*_ge)[_qp] * _bebar[_qp].deviatoric();
  ADReal s_trial_norm = s_trial.doubleContraction(s_trial);
  if (MooseUtils::absoluteFuzzyEqual(s_trial_norm, 0))
    s_trial_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
  s_trial_norm = sqrt(s_trial_norm);

  _Np[_qp] = sqrt(1.5) * s_trial / s_trial_norm;

  // Return mapping
  ADReal phi = computeResidual(s_trial_norm, delta_ep);
  if (phi > 0)
  {
    returnMappingSolve(s_trial_norm, delta_ep, _console);

    _ep[_qp] = _ep_old[_qp] + delta_ep;
    ADRankTwoTensor s = s_trial - (2.0 / 3.0) * (*_ge)[_qp] * (*_G)[_qp] * delta_ep *
                                      _bebar[_qp].trace() * _Np[_qp];

    ADReal J = _F[_qp].det();
    ADReal p = 0.5 * (*_K)[_qp] * (J * J - 1);
    ADRankTwoTensor tau = J >= 1.0 ? (*_ge)[_qp] * p * I2 + s : p * I2 + s;
    stress = tau / J;

    computeCorrectionTerm(s / (*_ge)[_qp] / (*_G)[_qp]);
  }
  else
  {
    _ep[_qp] = _ep_old[_qp];

    ADReal J = _F[_qp].det();
    ADReal p = 0.5 * (*_K)[_qp] * (J * J - 1);
    ADRankTwoTensor tau = J >= 1.0 ? (*_ge)[_qp] * p * I2 + s_trial : p * I2 + s_trial;
    stress = tau / J;
  }

  computeStrainEnergyDensity();
  _hardening_model->plasticEnergy(_ep[_qp]);
  _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 0);

  if (_t_step > 0)
  {
    _heat[_qp] = _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 1) * delta_ep / _dt;
    _heat[_qp] += _hardening_model->thermalConjugate(_ep[_qp]) * delta_ep / _dt;
  }
  else
    _heat[_qp] = 0;

  // psip_triax_raw: irreversible triaxiality-weighted plastic energy (no threshold).
  // Uses _triaxfunc from the previous NL iteration (lagged). On the first call it is
  // 0 so new_value is gated by the psi_a >= 0.1 guard.
  {
    const Real psi_a = MetaPhysicL::raw_value(_psip_active_ref[_qp]);
    const Real tf = MetaPhysicL::raw_value(_triaxfunc[_qp]);
    const Real old_raw = _psip_triax_raw_old[_qp];
    Real new_value = 0.0;
    if (psi_a >= 0.1 && tf > 0.0)
      new_value = psi_a / tf;
    if (new_value >= old_raw && new_value > 0.0)
      _psip_triax_raw[_qp] = _psip_active_ref[_qp] / _triaxfunc[_qp];
    else
      _psip_triax_raw[_qp] = ADReal(old_raw);
  }

  // psip_triax: threshold-gated version of psip_triax_raw.
  {
    const Real raw_value = MetaPhysicL::raw_value(_psip_triax_raw[_qp]);
    _psip_triax[_qp] = raw_value >= _psip_triax_threshold ? _psip_triax_raw[_qp] : ADReal(0.0);
  }

  // Compute stress triaxialities (Kirchhoff and Cauchy).
  {
    using std::sqrt;
    ADReal J = _F[_qp].det();
    ADRankTwoTensor tau = J * stress;

    ADReal tau_mean = tau.trace() / 3.0;
    ADRankTwoTensor s_tau = tau.deviatoric();
    ADReal s_tau_norm_sq = s_tau.doubleContraction(s_tau);
    if (MooseUtils::absoluteFuzzyEqual(s_tau_norm_sq, 0))
      s_tau_norm_sq.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
    _triaxiality_kirchhoff[_qp] = tau_mean / sqrt(1.5 * s_tau_norm_sq);

    ADReal sigma_mean = stress.trace() / 3.0;
    ADRankTwoTensor s_sigma = stress.deviatoric();
    ADReal s_sigma_norm_sq = s_sigma.doubleContraction(s_sigma);
    if (MooseUtils::absoluteFuzzyEqual(s_sigma_norm_sq, 0))
      s_sigma_norm_sq.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
    _triaxiality_cauchy[_qp] = sigma_mean / sqrt(1.5 * s_sigma_norm_sq);
  }

  // Gaussian triaxiality function: f(η) = 1 + (peak-1)*exp(-η²/(2σ²))
  // Peak at η=0, asymptotes to 1 at large |η|. σ² = 0.1803 gives σ ≈ 0.425.
  {
    ADReal eta = _triaxiality_cauchy[_qp];
    constexpr Real sigma_sq = 0.1803;
    _triaxfunc[_qp] = 1.0 + (_triax_gaussian_peak - 1.0) *
                                std::exp(MetaPhysicL::raw_value(-eta * eta / (2.0 * sigma_sq)));
  }

  _bebar_det[_qp] = _bebar[_qp].det();
}

Real
LargeDeformationJ2PlasticityBeBar::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                            const ADReal & delta_ep)
{
  return MetaPhysicL::raw_value(effective_trial_stress - std::sqrt(2.0 / 3.0) * (*_ge)[_qp] *
                                                             (*_G)[_qp] * delta_ep *
                                                             _bebar[_qp].trace());
}

ADReal
LargeDeformationJ2PlasticityBeBar::computeResidual(const ADReal & effective_trial_stress,
                                                   const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;
  return effective_trial_stress -
         std::sqrt(2.0 / 3.0) * (_hardening_model->plasticEnergy(ep, 1) +
                                 _hardening_model->plasticDissipation(delta_ep, ep, 1)) -
         std::sqrt(2.0 / 3.0) * (*_ge)[_qp] * (*_G)[_qp] * delta_ep * _bebar[_qp].trace();
}

ADReal
LargeDeformationJ2PlasticityBeBar::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                     const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;
  return -std::sqrt(2.0 / 3.0) * (_hardening_model->plasticEnergy(ep, 2) +
                                  _hardening_model->plasticDissipation(delta_ep, ep, 2)) -
         std::sqrt(2.0 / 3.0) * (*_ge)[_qp] * (*_G)[_qp] * _bebar[_qp].trace();
}

void
LargeDeformationJ2PlasticityBeBar::computeCorrectionTerm(const ADRankTwoTensor & devbebar)
{
  using std::cbrt;
  using std::sqrt;

  ADRankTwoTensor I2(RankTwoTensorTempl<ADReal>::initIdentity);

  ADReal a = devbebar(0, 0);
  ADReal b = devbebar(1, 1);
  ADReal c = devbebar(2, 2);
  ADReal d = devbebar(1, 2);
  ADReal e = devbebar(0, 2);
  ADReal h = devbebar(0, 1);

  ADReal A = a + b + c;
  ADReal B = a * b + a * c + b * c - d * d - e * e - h * h;
  ADReal C = a * b * c + 2.0 * d * e * h - a * d * d - b * e * e - c * h * h - 1.0;

  ADReal D = cbrt(
      -2 * A * A * A +
      3 * sqrt(3) *
          sqrt(4 * A * A * A * C - A * A * B * B - 18 * A * B * C + 4 * B * B * B + 27 * C * C) +
      9 * A * B - 27 * C);

  ADReal Ie_bar = D / 3 / cbrt(2) - cbrt(2) * (3 * B - A * A) / 3 / D - A / 3;
  _bebar[_qp] = devbebar + Ie_bar * I2;
}

void
LargeDeformationJ2PlasticityBeBar::computeStrainEnergyDensity()
{
  using std::log;
  ADReal J = _F[_qp].det();
  ADReal U = 0.5 * (*_K)[_qp] * (0.5 * (J * J - 1) - log(J));
  ADReal W = 0.5 * (*_G)[_qp] * (_bebar[_qp].trace() - 3.0);

  _psie_unsplit[_qp] = U + W;

  ADReal E_el_pos = J >= 1.0 ? U + W : W;
  ADReal E_el_neg = J >= 1.0 ? 0.0 : U;

  _psie_active_corr[_qp] = _apply_strain_energy_split ? E_el_pos : _psie_unsplit[_qp];
  _psie_corr[_qp] = (*_ge)[_qp] * E_el_pos + E_el_neg;
  _dpsie_dd_corr[_qp] = (*_dge_dd)[_qp] * _psie_active_corr[_qp];

  // Overwrite CNH's energy properties so that psie_active is correct for fracture driving.
  // CNH computes these from Fe, which BeBar does not update; the bebar-based values are correct.
  (*_psie_cnh)[_qp] = _psie_corr[_qp];
  (*_psie_active_cnh)[_qp] = _psie_active_corr[_qp];
  (*_dpsie_dd_cnh)[_qp] = _dpsie_dd_corr[_qp];
}
