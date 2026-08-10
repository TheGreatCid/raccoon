//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "LargeDeformationElasticityModel.h"
#include "DerivativeMaterialPropertyNameInterface.h"

/**
 * Isotropic compressible Neo-Hookean hyperelasticity with an alternative volumetric response.
 *
 * The volumetric part of the Mandel/Kirchhoff stress is
 *   tau_vol = 2 K (J - 1/J) (J + 1/J - 1) I,
 * which derives from the volumetric strain energy
 *   U(J) = 2 K (0.5 J^2 + 0.5 J^-2 - J - 1/J + 1)   (U(1) = 0, tau_vol = J dU/dJ I).
 * The isochoric (deviatoric) response is the standard Neo-Hookean G*dev(b) with the isochoric
 * distortional energy 0.5 G (tr(bbar) - 3), identical to CNHIsotropicElasticity.
 *
 * This model does not support any energy decomposition (no volumetric/deviatoric split) and does
 * not include the inversion barrier of CNHIsotropicElasticity.
 */
class CNHJay : public LargeDeformationElasticityModel,
               public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();

  CNHJay(const InputParameters & parameters);

  virtual ADRankTwoTensor computeMandelStress(const ADRankTwoTensor & Fe,
                                              const bool plasticity_update = false) override;

  const ADMaterialProperty<Real> & getK() const { return _K; }
  const ADMaterialProperty<Real> & getG() const { return _G; }
  const ADMaterialProperty<Real> & getDegradation() const { return _g; }
  const ADMaterialProperty<Real> & getDegradationDerivative() const { return _dg_dd; }

  ADMaterialProperty<Real> & getPsie() { return _psie; }
  ADMaterialProperty<Real> & getPsieActive() { return _psie_active; }
  ADMaterialProperty<Real> & getDpsieDD() { return _dpsie_dd; }

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
