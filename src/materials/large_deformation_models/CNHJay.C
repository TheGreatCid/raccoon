//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "CNHJay.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", CNHJay);

InputParameters
CNHJay::validParams()
{
  InputParameters params = LargeDeformationElasticityModel::validParams();
  params.addClassDescription(
      "Isotropic compressible Neo-Hookean hyperelasticity with the volumetric stress "
      "2*K*(J - 1/J)*(J + 1/J - 1)*I. No energy decomposition (split) is supported.");

  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");

  return params;
}

CNHJay::CNHJay(const InputParameters & parameters)
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
    _dg_dd(getADMaterialProperty<Real>(derivativePropertyName(_g_name, {_d_name})))
{
}

ADReal
CNHJay::volumetricKirchhoffPressure(const ADReal & J) const
{
  const ADReal Jinv = 1.0 / J;
  // Volumetric Kirchhoff pressure p (tau_vol = p I), = J dU/dJ for the m=1 Seth-Hill volumetric
  // energy below: p = (K/2) (J - 1/J) (J + 1/J - 1).
  return 0.5 * _K[_qp] * (J - Jinv) * (J + Jinv - 1.0);
}

ADReal
CNHJay::volumetricEnergy(const ADReal & J) const
{
  const ADReal Jinv = 1.0 / J;
  // Seth-Hill volumetric energy at m=1 (Garanger et al. 2026, Eq. 4): U = (K/4)[(J-1)^2+(1/J-1)^2].
  // U(1) = 0, and tau_vol = J dU/dJ I (see volumetricKirchhoffPressure).
  return 0.25 * _K[_qp] * ((J - 1.0) * (J - 1.0) + (Jinv - 1.0) * (Jinv - 1.0));
}

ADRankTwoTensor
CNHJay::computeMandelStress(const ADRankTwoTensor & Fe, const bool plasticity_update)
{
  using std::pow;
  using std::sqrt;

  // We use the left Cauchy-Green strain
  ADRankTwoTensor strain;
  if (plasticity_update)
  {
    ADRankTwoTensor expFe = RaccoonUtils::exp(Fe);
    strain = expFe * expFe.transpose();
  }
  else
    strain = Fe * Fe.transpose();

  const ADReal J = sqrt(strain.det());

  const ADRankTwoTensor I2(ADRankTwoTensor::initIdentity);

  // Volumetric Mandel/Kirchhoff stress p*I plus the standard Neo-Hookean deviatoric G*dev(b). The
  // volumetric term is sourced from volumetricKirchhoffPressure so the elastic path and any
  // plasticity model that queries this model share a single definition.
  ADRankTwoTensor stress_intact =
      volumetricKirchhoffPressure(J) * I2 + _G[_qp] * strain.deviatoric();
  ADRankTwoTensor stress = _g[_qp] * stress_intact;

  // If plasticity_update == false, then we are not in the middle of a plasticity update, hence we
  // compute the strain energy density.
  if (!plasticity_update)
  {
    ADRankTwoTensor strain_bar = pow(J, -2. / 3.) * strain;
    ADReal U = volumetricEnergy(J);
    ADReal W = 0.5 * _G[_qp] * (strain_bar.trace() - 3.0);
    _psie_active[_qp] = U + W;
    _psie[_qp] = _g[_qp] * _psie_active[_qp];
    _dpsie_dd[_qp] = _dg_dd[_qp] * _psie_active[_qp];
  }

  return stress;
}
