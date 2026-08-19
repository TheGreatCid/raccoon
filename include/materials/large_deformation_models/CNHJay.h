//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "LargeDeformationElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"
#include "CNHElasticityInterface.h"

/**
 * Isotropic compressible Neo-Hookean hyperelasticity with an alternative volumetric response.
 *
 * The volumetric response is the Seth-Hill volumetric term at m=1 (Garanger et al. 2026, Eq. 4):
 *   U(J)     = (K/4) [ (J - 1)^2 + (1/J - 1)^2 ]        (U(1) = 0)
 *   tau_vol  = J dU/dJ I = (K/2) (J - 1/J) (J + 1/J - 1) I.
 * The isochoric (deviatoric) response is the standard Neo-Hookean G*dev(b) with the isochoric
 * distortional energy 0.5 G (tr(bbar) - 3), identical to CNHIsotropicElasticity (this is NOT the
 * Seth-Hill deviatoric term of Eq. 4).
 *
 * This model does not support any energy decomposition (no volumetric/deviatoric split) and does
 * not include the inversion barrier of CNHIsotropicElasticity.
 */
class CNHJay : public LargeDeformationElasticityModel,
               public DerivativeMaterialPropertyNameInterface,
               public CNHElasticityInterface
{
public:
  static InputParameters validParams();

  CNHJay(const InputParameters & parameters);

  virtual ADRankTwoTensor computeMandelStress(const ADRankTwoTensor & Fe,
                                              const bool plasticity_update = false) override;

  const ADMaterialProperty<Real> & getK() const override { return _K; }
  const ADMaterialProperty<Real> & getG() const override { return _G; }
  const ADMaterialProperty<Real> & getDegradation() const override { return _g; }
  const ADMaterialProperty<Real> & getDegradationDerivative() const override { return _dg_dd; }

  ADMaterialProperty<Real> & getPsie() override { return _psie; }
  ADMaterialProperty<Real> & getPsieActive() override { return _psie_active; }
  ADMaterialProperty<Real> & getDpsieDD() override { return _dpsie_dd; }

  ADReal volumetricKirchhoffPressure(const ADReal & J) const override;
  ADReal volumetricEnergy(const ADReal & J) const override;
  bool supportsEnergySplit() const override { return false; }

protected:
  /// The bulk modulus
  const ADMaterialProperty<Real> & _K;

  /// The shear modulus
  const ADMaterialProperty<Real> & _G;

  /// Name of the phase-field variable
  const VariableName _d_name;

  // @{ Strain energy density and its derivative w/r/t damage
  const MaterialPropertyName _psie_name;
  ADMaterialProperty<Real> & _psie;
  ADMaterialProperty<Real> & _psie_active;
  ADMaterialProperty<Real> & _dpsie_dd;
  // @}

  // @{ The degradation function and its derivative w/r/t damage
  const MaterialPropertyName _g_name;
  const ADMaterialProperty<Real> & _g;
  const ADMaterialProperty<Real> & _dg_dd;
  // @}
};
