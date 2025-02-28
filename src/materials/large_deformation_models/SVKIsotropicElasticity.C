//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SVKIsotropicElasticity.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", SVKIsotropicElasticity);

InputParameters
SVKIsotropicElasticity::validParams()
{
  InputParameters params = LargeDeformationElasticityModel::validParams();
  params.addClassDescription("Isotropic Compressible Neo-Hookean hyperelasticity model.");

  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");
  params.addParam<MooseEnum>(
      "decomposition", MooseEnum("NONE", "NONE"), "The decomposition method");

  return params;
}

SVKIsotropicElasticity::SVKIsotropicElasticity(const InputParameters & parameters)
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
SVKIsotropicElasticity::computeMandelStress(const ADRankTwoTensor & Fe,
                                            const bool plasticity_update)
{
  ADRankTwoTensor stress;

  if (_decomposition == Decomposition::none)
    stress = computeMandelStressNoDecomposition(Fe, plasticity_update);
  else
    paramError("decomposition", "Unsupported decomposition type.");

  return stress;
}

ADRankTwoTensor
SVKIsotropicElasticity::computeMandelStressNoDecomposition(const ADRankTwoTensor & Fe,
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

  //   ADReal J = std::sqrt(strain.det());
  auto s2 = strain * strain;
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);
  auto lambda = (_K[_qp] - 2 / 3 * _G[_qp]);
  ADRankTwoTensor stress_intact =
      _G[_qp] * (s2 - strain) + 0.5 * lambda * (strain.trace() - 3) * strain;
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  // If plasticity_update == false, then we are not in the middle of a plasticity update, hence we
  // compute the strain energy density
  if (!plasticity_update)
  {

    ADReal U = 0.5 * lambda * strain.trace() * strain.trace() + _G[_qp] * (s2).trace();
    _psie_active[_qp] = U;
    _psie[_qp] = _g[_qp] * _psie_active[_qp];
    _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];
  }

  return stress;
}
