//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "EigenADReal.h"
#include "LargeDeformationJ2PlasticityCorrection.h"
#include "LargeDeformationPlasticityModel.h"
#include "MooseError.h"
#include "MooseTypes.h"
#include "RaccoonUtils.h"
#include "RankTwoTensorForward.h"
#include <cmath>

registerMooseObject("raccoonApp", LargeDeformationJ2PlasticityCorrection);

InputParameters
LargeDeformationJ2PlasticityCorrection::validParams()
{
  InputParameters params = LargeDeformationPlasticityModel::validParams();

  params.addClassDescription("Large deformation $J_2$ plasticity. The exponential constitutive "
                             "update is used to update the plastic deformation.");
  params.addRequiredParam<MaterialPropertyName>("K", "K");
  params.addRequiredParam<MaterialPropertyName>("G", "G");
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "elastic_degradation_function", "g", "The elastic degradation function");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density_corr",
      "psie_corr",
      "Name of the strain energy density computed by this material model");

  return params;
}

LargeDeformationJ2PlasticityCorrection::LargeDeformationJ2PlasticityCorrection(
    const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _bebar(declareADProperty<RankTwoTensor>(prependBaseName("be_bar"))),
    _bebar_old(getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("be_bar"))),
    // _f(declareADProperty<RankTwoTensor>(prependBaseName("incremental_deformation_gradient"))),
    _G(getADMaterialPropertyByName<Real>("G")),
    _K(getADMaterialPropertyByName<Real>("K")),
    _F_old(getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _F(getADMaterialPropertyByName<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _d_name(getVar("phase_field", 0)->name()),
    // The degradation function and its derivatives
    _ge_name(prependBaseName("elastic_degradation_function", true)),
    _ge(getADMaterialProperty<Real>(_ge_name)),
    _dge_dd(getADMaterialProperty<Real>(derivativePropertyName(_ge_name, {_d_name}))),

    // The strain energy density and its derivatives
    _psie_name(prependBaseName("strain_energy_density_corr", true)),
    _psie_corr(declareADProperty<Real>(_psie_name)),
    _psie_active_corr(declareADProperty<Real>(_psie_name + "_active")),
    _dpsie_dd_corr(declareADProperty<Real>(derivativePropertyName(_psie_name, {_d_name})))

{
  _check_range = true;
}

void
LargeDeformationJ2PlasticityCorrection::initQpStatefulProperties()
{
  _Fp[_qp].setToIdentity();
  _ep[_qp] = 0;
  _bebar[_qp].setToIdentity();
  std::vector<std::string> indices = {"x", "y", "z"};
  if (_recover)
  {
    _ep[_qp] =
        _solution_object_ptr->pointValue(_t,
                                         _current_elem->true_centroid(),
                                         "effective_plastic_strain_" + std::to_string(_qp + 1),
                                         nullptr);

    for (int i_ind = 0; i_ind < 3; i_ind++)
      for (int j_ind = 0; j_ind < 3; j_ind++)
      {
        _bebar[_qp](i_ind, j_ind) = _solution_object_ptr->pointValue(
            _t,
            _current_elem->true_centroid(),
            "be_bar_" + indices[i_ind] + indices[j_ind] + "_" + std::to_string(_qp + 1),
            nullptr);
      }
  }
};

void
LargeDeformationJ2PlasticityCorrection::updateState(ADRankTwoTensor & stress,
                                                    ADRankTwoTensor & /*Fe*/)
{
  ADRankTwoTensor I2;
  I2.setToIdentity();

  // Assume no plastic increment
  ADReal delta_ep = 0;

  // Obtain incremental f
  ADRankTwoTensor f = _F[_qp] * _F_old[_qp].inverse();

  // Compute fbar
  ADRankTwoTensor fbar = f / std::cbrt(f.det());

  // Not sure if this actually helps
  if (_t_step > 0 || _recover == false)
    // Compute bebar_trial
    _bebar[_qp] = fbar * _bebar_old[_qp] * fbar.transpose();

  // Compute trial stress
  // Will need to do a better implementation of this
  ADRankTwoTensor s_trial = _G[_qp] * _ge[_qp] * _bebar[_qp].deviatoric();
  ADReal s_trial_norm = s_trial.doubleContraction(s_trial);
  if (MooseUtils::absoluteFuzzyEqual(s_trial_norm, 0))
    s_trial_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;

  // Norm of trial stress
  s_trial_norm = std::sqrt(s_trial_norm);

  // Check for plastic loading
  _Np[_qp] = std::sqrt(1.5) * s_trial / s_trial_norm;

  // Return mapping
  ADReal phi = computeResidual(s_trial_norm, delta_ep);
  if (phi > 0)
  {
    returnMappingSolve(s_trial_norm, delta_ep, _console);

    // Update stress
    ADRankTwoTensor s = s_trial - _ge[_qp] * _G[_qp] * delta_ep * _bebar[_qp].trace() * _Np[_qp];
    _ep[_qp] = _ep_old[_qp] + delta_ep;

    // Updating Kirchoff stress
    ADReal J = _F[_qp].det();
    ADReal p = 0.5 * _K[_qp] * (J * J - 1);
    stress = J >= 1.0 ? _ge[_qp] * p * I2 + s : p * I2 + s;

    ADRankTwoTensor devbebar = s / _ge[_qp] / _G[_qp];
    computeCorrectionTerm(devbebar);
  }
  else
  {
    // Updating Kirchoff stress
    ADReal J = _F[_qp].det();
    ADReal p = 0.5 * _K[_qp] * (J * J - 1);
    stress = J * p * I2 + s_trial;
    _ep[_qp] = _ep_old[_qp];
  }
  computeStrainEnergyDensity();
  _hardening_model->plasticEnergy(_ep[_qp]);
  _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 0);
  // if (_current_elem->id() == 1)
  // {
  //   std::cout << "=================" << std::endl;
  //   MetaPhysicL::raw_value(_bebar_old[_qp]).print();

  //   MetaPhysicL::raw_value(_bebar[_qp]).print();
  //   std::cout << MetaPhysicL::raw_value(delta_ep) << std::endl;
  //   std::cout << "=================" << std::endl;
  // }
}

Real
LargeDeformationJ2PlasticityCorrection::computeReferenceResidual(
    const ADReal & effective_trial_stress, const ADReal & delta_ep)
{
  return MetaPhysicL::raw_value(effective_trial_stress -
                                std::sqrt(3.0 / 2.0) * _G[_qp] * delta_ep * _bebar[_qp].trace());
}

ADReal
LargeDeformationJ2PlasticityCorrection::computeResidual(const ADReal & effective_trial_stress,
                                                        const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;
  return effective_trial_stress - std::sqrt(2.0 / 3.0) * _hardening_model->plasticEnergy(ep, 1) -
         std::sqrt(3.0 / 2.0) * _ge[_qp] * _G[_qp] * delta_ep * _bebar[_qp].trace();
}

ADReal
LargeDeformationJ2PlasticityCorrection::computeDerivative(const ADReal & /*effective_trial_stress*/
                                                          ,
                                                          const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;
  return -std::sqrt(2.0 / 3.0) * _hardening_model->plasticEnergy(ep, 2) -
         std::sqrt(3.0 / 2.0) * _ge[_qp] * _G[_qp] * _bebar[_qp].trace();
}

void
LargeDeformationJ2PlasticityCorrection::computeCorrectionTerm(const ADRankTwoTensor & devbebar)
{

  // Identity tensor
  ADRankTwoTensor I2(RankTwoTensorTempl<ADReal>::initIdentity);

  Real a = MetaPhysicL::raw_value(devbebar(0, 0));
  Real b = MetaPhysicL::raw_value(devbebar(1, 1));
  Real c = MetaPhysicL::raw_value(devbebar(2, 2));
  Real d = MetaPhysicL::raw_value(devbebar(1, 2));
  Real e = MetaPhysicL::raw_value(devbebar(0, 2));
  Real h = MetaPhysicL::raw_value(devbebar(0, 1));

  Real A = a + b + c;
  Real B = a * b + a * c + b * c - d * d - e * e - h * h;
  Real C = a * b * c + 2.0 * d * e * h - a * d * d - b * e * e - c * h * h - 1.0;

  Real D = std::cbrt(-2 * A * A * A +
                     3 * std::sqrt(3) *
                         std::sqrt(4 * A * A * A * C - A * A * B * B - 18 * A * B * C +
                                   4 * B * B * B + 27 * C * C) +
                     9 * A * B - 27 * C);

  ADReal Ie_bar = D / 3 / std::cbrt(2) - std::cbrt(2) * (3 * B - A * A) / 3 / D - A / 3;
  _bebar[_qp] = devbebar + Ie_bar * I2;
}

void
LargeDeformationJ2PlasticityCorrection::computeStrainEnergyDensity()
{
  ADReal J = _F[_qp].det();
  ADReal U = 0.5 * _K[_qp] * (0.5 * (J * J - 1) - std::log(J));
  ADReal W = 0.5 * _G[_qp] * (_bebar[_qp].trace() - 3.0);

  ADReal E_el_pos = J >= 1.0 ? U + W : W;
  ADReal E_el_neg = J >= 1.0 ? 0.0 : U;

  _psie_active_corr[_qp] = E_el_pos;
  _psie_corr[_qp] = _ge[_qp] * E_el_pos + E_el_neg;

  _dpsie_dd_corr[_qp] = _dge_dd[_qp] * _psie_active_corr[_qp];
}