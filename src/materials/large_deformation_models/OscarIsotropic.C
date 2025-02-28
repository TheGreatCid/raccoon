//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "EigenADReal.h"
#include "OscarIsotropic.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", OscarIsotropic);

InputParameters
OscarIsotropic::validParams()
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
  params.addParam<double>("alpha", 1 / 5, "alpha term");
  return params;
}

OscarIsotropic::OscarIsotropic(const InputParameters & parameters)
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

    _decomposition(getParam<MooseEnum>("decomposition").getEnum<Decomposition>()),

    _alpha(getParam<double>("alpha"))
{
}

ADRankTwoTensor
OscarIsotropic::computeMandelStress(const ADRankTwoTensor & Fe, const bool plasticity_update)
{
  ADRankTwoTensor stress;

  if (_decomposition == Decomposition::none)
    stress = computeMandelStressNoDecomposition(Fe, plasticity_update);
  else
    paramError("decomposition", "Unsupported decomposition type.");

  return stress;
}

ADRankTwoTensor
OscarIsotropic::computeMandelStressNoDecomposition(const ADRankTwoTensor & Fe,
                                                   const bool plasticity_update)
{
  // We use the left Cauchy-Green strain
  ADRankTwoTensor strain;
  if (plasticity_update)
  {
    ADRankTwoTensor expFe = RaccoonUtils::exp(Fe);
    strain = expFe * expFe.transpose();
    if (_t == 2040)
    {
      MetaPhysicL::raw_value(strain).print();
    }
  }
  else
    strain = Fe * Fe.transpose();

  ADReal J = std::sqrt(strain.det());
  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  ADRankTwoTensor stress_intact =
      -std::pow(3, 1 - _alpha) * _G[_qp] * std::pow((strain.inverse()).trace(), _alpha - 1) *
          strain.inverse() +
      ((_K[_qp] - 2 / 3 * _alpha * _G[_qp]) * J * (J - 1) + _G[_qp]) * I2;

  // ADRankTwoTensor stress_intact =
  //     std::pow(3, 1 - _alpha) * _G[_qp] * std::pow(strain.trace(), _alpha - 1) * strain +
  //     ((_K[_qp] - 2 / 3 * _alpha * _G[_qp]) * J * (J - 1) - _G[_qp]) * I2;
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  // If plasticity_update == false, then we are not in the middle of a plasticity update, hence we
  // compute the strain energy density
  if (!plasticity_update)
  {

    ADReal U = 1;
    _psie_active[_qp] = U;
    _psie[_qp] = _g[_qp] * _psie_active[_qp];
    _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];
  }

  return stress;
}
