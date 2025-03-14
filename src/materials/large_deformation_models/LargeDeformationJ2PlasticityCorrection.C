//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "EigenADReal.h"
#include "LargeDeformationJ2PlasticityCorrection.h"
#include "LargeDeformationPlasticityModel.h"
#include "MooseTypes.h"
#include "RaccoonUtils.h"
#include "RankTwoTensorForward.h"

registerMooseObject("raccoonApp", LargeDeformationJ2PlasticityCorrection);

InputParameters
LargeDeformationJ2PlasticityCorrection::validParams()
{
  InputParameters params = LargeDeformationPlasticityModel::validParams();

  params.addClassDescription("Large deformation $J_2$ plasticity. The exponential constitutive "
                             "update is used to update the plastic deformation.");
  params.addRequiredParam<MaterialPropertyName>("K", "K");
  params.addRequiredParam<MaterialPropertyName>("G", "G");

  return params;
}

LargeDeformationJ2PlasticityCorrection::LargeDeformationJ2PlasticityCorrection(
    const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    _bebar(declareADProperty<RankTwoTensor>(prependBaseName("be_bar"))),
    _bebar_old(getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("be_bar"))),
    // _f(declareADProperty<RankTwoTensor>(prependBaseName("incremental_deformation_gradient"))),
    _G(getADMaterialPropertyByName<Real>("G")),
    _K(getADMaterialPropertyByName<Real>("K")),
    _F_old(getMaterialPropertyOldByName<RankTwoTensor>(prependBaseName("deformation_gradient"))),
    _F(getADMaterialPropertyByName<RankTwoTensor>(prependBaseName("deformation_gradient")))
{
  _check_range = true;
}

void
LargeDeformationJ2PlasticityCorrection::initQpStatefulProperties()
{
  _Fp[_qp].setToIdentity();
  _ep[_qp] = 0;
  _bebar[_qp].setToIdentity();
};

void
LargeDeformationJ2PlasticityCorrection::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
{
  ADRankTwoTensor I2;
  I2.setToIdentity();
  // Assume no plastic increment
  ADReal delta_ep = 0;

  // Obtain incremental f
  ADRankTwoTensor f = _F[_qp] * _F_old[_qp].inverse();

  // Compute fbar
  ADRankTwoTensor fbar = std::pow(f.det(), -1.0 / 3.0) * f;

  // Compute bebar_trial
  _bebar[_qp] = fbar * _bebar_old[_qp] * fbar.transpose();

  // Compute trial stress
  // Will need to do a better implementation of this
  ADRankTwoTensor s_trial = _G[_qp] * _bebar[_qp].deviatoric();
  ADReal s_trial_norm = s_trial.doubleContraction(s_trial);
  if (MooseUtils::absoluteFuzzyEqual(s_trial_norm, 0))
    s_trial_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;

  // Norm of trial stress
  s_trial_norm = std::sqrt(1.5 * s_trial_norm);

  // Check for plastic loading
  _Np[_qp] = 1.5 * s_trial / s_trial_norm;
  // Return mapping
  ADReal phi = computeResidual(s_trial_norm, delta_ep);
  if (phi > 0)
    returnMappingSolve(s_trial_norm, delta_ep, _console);

  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];
  _ep[_qp] = _ep_old[_qp] + delta_ep;

  // Update stress and energy
  Fe = Fe * delta_Fp.inverse();

  // Update stress
  ADReal Iebar = 1.0 / 3.0 * _bebar[_qp].trace();
  ADReal mubar = _G[_qp] * Iebar;
  ADRankTwoTensor s = s_trial - 2 * mubar * delta_ep * _Np[_qp];

  ADReal J = Fe.det();

  // Updating Kirchoff stress
  //   ADReal p = _K[_qp] / 2 * (J - 1 / J);
  ADReal p = 0.5 * _K[_qp] * (J * J - 1);
  stress = J * p * I2 + s;
  ADRankTwoTensor devbebar = stress.deviatoric() / _G[_qp];
  computeCorrectionTerm(devbebar);
}

Real
LargeDeformationJ2PlasticityCorrection::computeReferenceResidual(
    const ADReal & effective_trial_stress, const ADReal & delta_ep)
{
  ADReal Iebar = 1.0 / 3.0 * _bebar[_qp].trace();
  ADReal mubar = _G[_qp] * Iebar;
  return MetaPhysicL::raw_value(effective_trial_stress - 2 * mubar * delta_ep);
}

ADReal
LargeDeformationJ2PlasticityCorrection::computeResidual(const ADReal & effective_trial_stress,
                                                        const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;

  ADReal Iebar = 1.0 / 3.0 * _bebar[_qp].trace();
  ADReal mubar = _G[_qp] * Iebar;
  return effective_trial_stress - std::sqrt(2.0 / 3.0) * _hardening_model->plasticEnergy(ep, 1) -
         2 * mubar * delta_ep;
}

ADReal
LargeDeformationJ2PlasticityCorrection::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                          const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;
  ADReal Iebar = 1 / 3 * _bebar[_qp].trace();
  ADReal mubar = _G[_qp] * Iebar;
  return -std::sqrt(2.0 / 3.0) * _hardening_model->plasticEnergy(ep, 2) - 2 * mubar;
}

void
LargeDeformationJ2PlasticityCorrection::computeCorrectionTerm(const ADRankTwoTensor & devbebar)
{

  std::vector<ADReal> d;
  ADRankTwoTensor V;
  devbebar.symmetricEigenvaluesEigenvectors(d, V);

  // NR its: solving for det(bebar + alpha) = 1
  ADReal alpha = 1;
  ADReal i = 0;
  ADReal res = 1;
  ADReal dres;
  while (std::abs(res) > 1e-10 && i <= 50)
  {
    res = computeCorrectionResidual(d, alpha);
    dres = computeCorrectionDerivative(d, alpha);

    alpha = alpha - res / dres;

    res = computeCorrectionResidual(d, alpha);
    i += 1;
  }

  // Reconstructing bebar
  ADRankTwoTensor d_diag;
  d_diag.setToIdentity();
  d_diag(0, 0) = d[0] + alpha;
  d_diag(1, 1) = d[1] + alpha;
  d_diag(2, 2) = d[2] + alpha;
  _bebar[_qp] = V * d_diag * V.transpose();
}

ADReal
LargeDeformationJ2PlasticityCorrection::computeCorrectionResidual(const std::vector<ADReal> d,
                                                                  const ADReal & alpha)
{
  return (d[0] + alpha) * (d[1] + alpha) * (d[2] + alpha) - 1;
}

ADReal
LargeDeformationJ2PlasticityCorrection::computeCorrectionDerivative(const std::vector<ADReal> d,
                                                                    const ADReal & alpha)
{
  return (d[1] + alpha) * (d[2] + alpha) + (d[0] + alpha) * (d[2] + alpha) +
         (d[0] + alpha) * (d[1] + alpha);
}