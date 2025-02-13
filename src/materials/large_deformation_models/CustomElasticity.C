//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "CustomElasticity.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", CustomElasticity);

InputParameters
CustomElasticity::validParams()
{
  InputParameters params = LargeDeformationElasticityModel::validParams();
  params.addClassDescription("Isotropic Compressible Neo-Hookean elasticity model.");

  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");
  params.addParam<MooseEnum>("decomposition", MooseEnum("NONE"), "The decomposition method");

  return params;
}

CustomElasticity::CustomElasticity(const InputParameters & parameters)
  : LargeDeformationElasticityModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _K(getADMaterialPropertyByName<Real>(prependBaseName("bulk_modulus", true))),
    _G(getADMaterialPropertyByName<Real>(prependBaseName("shear_modulus", true))),

    _d_name(getVar("phase_field", 0)->name()),

    // The strain energy density and its derivatives
    _psie_name(prependBaseName("strain_energy_density", true)),
    _psie(declareADProperty<Real>(_psie_name)),
    _psie_active(declareADProperty<Real>(_psie_name + "_active")),
    _dpsie_dd(declareADProperty<Real>(derivativePropertyName(_psie_name, {_d_name}))),

    // The degradation function and its derivatives
    _g_name(prependBaseName("degradation_function", true)),
    _g(getADMaterialProperty<Real>(_g_name)),
    _dg_dd(getADMaterialProperty<Real>(derivativePropertyName(_g_name, {_d_name}))),

    _decomposition(getParam<MooseEnum>("decomposition").getEnum<Decomposition>())
{
}

ADRankTwoTensor
CustomElasticity::computeMandelStress(const ADRankTwoTensor & Fe, const bool plasticity_update)
{
  ADRankTwoTensor stress;

  if (_decomposition == Decomposition::none)
    stress = computeMandelStressNoDecomposition(Fe, plasticity_update);
  else
    paramError("decomposition", "Unsupported decomposition type.");

  return stress;
}

ADRankTwoTensor
CustomElasticity::computeMandelStressNoDecomposition(const ADRankTwoTensor & Fe,
                                                     const bool plasticity_update)
{
  // We use the left Cauchy-Green strain
  ADRankTwoTensor strain;
  if (plasticity_update)
  {
    ADRankTwoTensor expFe = RaccoonUtils::exp(Fe);
    strain = expFe * expFe.transpose();
  }
  else
    strain = Fe * Fe.transpose();

  ADReal J = std::sqrt(strain.det());

  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  // ADRankTwoTensor stress_intact =
  //     (_K[_qp] * (J - 1) * J + _G[_qp] * ((1 / 3) * strain.trace() - 1 - 2 / 3 * (J - 1) * J)) *
  //         I2 +
  //     _G[_qp] * strain.deviatoric();
  ADRankTwoTensor stress_intact =
      _G[_qp] * (strain - I2) + (_K[_qp] - (2.0 / 3.0) * _G[_qp]) * (J - 1) * J * I2;
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  // // Compute the first Piola-Kirchoff stress
  // ADRankTwoTensor P = computeFirstPKStressNoDecomposition(Fe, plasticity_update);

  // Compute the undegraded and degraded Mandel stress (here M = kirchoff stress)
  // ADRankTwoTensor stress_intact = P * Fe.transpose();
  // ADRankTwoTensor stress = _g[_qp] * stress_intact;

  if (!plasticity_update)
  {
    // Compute the strain energy density
    ADRankTwoTensor C = Fe.transpose() * Fe;
    ADReal J = std::sqrt(C.det());
    ADReal lambda = _K[_qp] - (2.0 / 3.0) * _G[_qp];
    ADReal W_vol = -_G[_qp] * std::log(J) + 0.5 * lambda * std::pow((J - 1), 2);
    ADReal I1 = C.trace(); // First invariant of C
    ADReal W_dev = 0.5 * _G[_qp] * (I1 - 3.0);
    _psie_active[_qp] = W_vol + W_dev;
    _psie[_qp] = _g[_qp] * _psie_active[_qp];
    _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];
  }

  return stress;
}

// ADRankTwoTensor
// CustomElasticity::computeFirstPKStressNoDecomposition(const ADRankTwoTensor & Fe,
//                                                       const bool plasticity_update)
// {

//   // We use the right Cauchy-Green strain

//   ADReal J = std::sqrt(strain.det());

//   // Compute the first Piola-Kirchoff Stress
//   ADRankTwoTensor Fe_invT = Fe.inverse().transpose();
//   ADReal lambda = _K[_qp] - (2.0 / 3.0) * _G[_qp];
//   ADRankTwoTensor P = _G[_qp] * (Fe - Fe_invT) + lambda * (J - 1) * J * Fe_invT;

//   return P;
// }