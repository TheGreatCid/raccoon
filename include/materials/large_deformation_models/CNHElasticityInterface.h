//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "MaterialProperty.h"

/**
 * Interface exposing the elastic properties and degraded strain-energy densities that the
 * bebar J2 plasticity update (LargeDeformationJ2PlasticityBeBar) needs from its companion
 * Neo-Hookean elasticity model. Both CNHIsotropicElasticity and CNHJay implement it, so the
 * plasticity model can source these quantities from either without depending on a concrete type.
 */
class CNHElasticityInterface
{
public:
  virtual ~CNHElasticityInterface() = default;

  virtual const ADMaterialProperty<Real> & getK() const = 0;
  virtual const ADMaterialProperty<Real> & getG() const = 0;
  virtual const ADMaterialProperty<Real> & getDegradation() const = 0;
  virtual const ADMaterialProperty<Real> & getDegradationDerivative() const = 0;

  virtual ADMaterialProperty<Real> & getPsie() = 0;
  virtual ADMaterialProperty<Real> & getPsieActive() = 0;
  virtual ADMaterialProperty<Real> & getDpsieDD() = 0;

  /**
   * The volumetric part of the Kirchhoff/Mandel stress is p*I, with p returned here. This lets a
   * plasticity model (e.g. LargeDeformationJ2PlasticityBeBar) build the volumetric response from the
   * companion elasticity model's own formulation instead of hard-coding one. Uses the model's bulk
   * modulus at the model's current quadrature point.
   */
  virtual ADReal volumetricKirchhoffPressure(const ADReal & J) const = 0;

  /// The volumetric strain energy density U(J), consistent with volumetricKirchhoffPressure (p = J dU/dJ).
  virtual ADReal volumetricEnergy(const ADReal & J) const = 0;

  /**
   * Whether this model's energy admits a volumetric/deviatoric (tension/compression) split. Models
   * that define a single unsplit response (e.g. CNHJay) return false, and a plasticity model should
   * then not apply an energy split regardless of its own setting.
   */
  virtual bool supportsEnergySplit() const = 0;

  /**
   * Optional un-degraded inversion-barrier Kirchhoff pressure that a plasticity model should add on
   * top of the volumetric stress in (nearly) collapsed damaged elements to resist inversion. Returns
   * 0 when the model has no barrier or the current point is not damaged enough. Default: no barrier.
   * See CNHIsotropicElasticity for the enabling parameters.
   */
  virtual ADReal inversionBarrierPressure(const ADReal & /*J*/) const { return 0.0; }

  /// Un-degraded strain energy density of the inversion barrier (0 if none), consistent with
  /// inversionBarrierPressure. It is added to the total energy only and does not drive fracture.
  virtual ADReal inversionBarrierEnergy(const ADReal & /*J*/) const { return 0.0; }
};
