//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "AuxKernel.h"

/**
 * Kinetic energy density in the "product form" that the leapfrog
 * (central-difference) time integrator conserves exactly:
 *
 *     KE^n = 1/2 rho ( v^{n+1/2} . v^{n-1/2} )
 *
 * The velocity vector produced by ExplicitMixedOrder is the leapfrog
 * half-step velocity v^{n+1/2}.  Mirrored into the coupled velocity
 * AuxVariables (via CoupledTimeDerivativeAux at TIMESTEP_END), its current
 * value is v^{n+1/2} and its OLD value (previous timestep) is v^{n-1/2}.
 *
 * Using the product of the two surrounding half-step velocities -- rather
 * than 1/2 rho |v^{n+1/2}|^2 -- yields a discrete total energy (KE + strain)
 * that stays flat in time, making it the right quantity for checking whether
 * the scheme conserves energy.  Contrast KineticEnergyAux, which squares a
 * single half-step velocity and is therefore offset dt/2 from the strain
 * energy.
 *
 * The density is read from a material property so the energy stays correct
 * for a spatially varying rho.
 */
class LeapfrogKineticEnergyAux : public AuxKernel
{
public:
  static InputParameters validParams();
  LeapfrogKineticEnergyAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Base name of the material system (prefix for the density property)
  const std::string _base_name;

  /// Density rho, from a material property (supports spatially varying density)
  const MaterialProperty<Real> & _density;

  /// Current half-step velocity components: v^{n+1/2}
  const VariableValue & _vel_x;
  const VariableValue & _vel_y;
  const VariableValue & _vel_z;

  /// Previous half-step velocity components: v^{n-1/2} (OLD state of vel_*)
  const VariableValue & _vel_x_old;
  const VariableValue & _vel_y_old;
  const VariableValue & _vel_z_old;
};
