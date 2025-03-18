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

  // Compute bebar_trial
  _bebar[_qp] = fbar * _bebar_old[_qp] * fbar.transpose();

  // Compute trial stress
  // Will need to do a better implementation of this
  ADRankTwoTensor s_trial = _G[_qp] * _bebar[_qp].deviatoric();
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
    returnMappingSolve(s_trial_norm, delta_ep, _console);

  // Update stress
  ADRankTwoTensor s = s_trial - _G[_qp] * delta_ep * _bebar[_qp].trace() * _Np[_qp];
  _ep[_qp] = _ep_old[_qp] + delta_ep * 2;

  // Updating Kirchoff stress
  ADReal J = _F[_qp].det();
  ADReal p = 0.5 * _K[_qp] * (J * J - 1);
  stress = J * p * I2 + s;
  ADRankTwoTensor devbebar = s / _G[_qp];
  computeCorrectionTerm(devbebar);
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
         std::sqrt(3.0 / 2.0) * _G[_qp] * delta_ep * _bebar[_qp].trace();
}

ADReal
LargeDeformationJ2PlasticityCorrection::computeDerivative(const ADReal & /*effective_trial_stress*/
                                                          ,
                                                          const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;
  return -std::sqrt(2.0 / 3.0) * _hardening_model->plasticEnergy(ep, 2) -
         std::sqrt(3.0 / 2.0) * _G[_qp] * _bebar[_qp].trace();
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
