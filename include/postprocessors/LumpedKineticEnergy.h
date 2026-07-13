//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "GeneralPostprocessor.h"

/**
 * Kinetic energy measured in the integrator's OWN lumped-mass norm:
 *
 *     KE = 1/2 sum_i m_i (u_dot_i)^2
 *
 * where m_i is the lumped mass diagonal that ExplicitMixedOrder(HRZ/Recover)
 * assembles into its "mass_matrix_lumped" system vector, and u_dot_i is the
 * integrator's half-step velocity (solutionUDot).  This is exactly the
 * kinetic energy the central-difference scheme conserves against -- unlike an
 * ElementIntegralVariablePostprocessor of 1/2 rho v^2, which uses consistent
 * (Gauss) quadrature and a user-chosen density.
 *
 * Diagnostic use: run it in both the reference and the recover legs.  Because
 * the recovered half-step velocity matches the reference to ~1e-4, any larger
 * gap in this quantity at the handoff isolates the difference in the lumped
 * mass M itself (the adj_density / det(F) reconstruction error).
 */
class LumpedKineticEnergy : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  LumpedKineticEnergy(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual PostprocessorValue getValue() const override;

protected:
  /// Nonlinear system that owns the lumped mass vector and solutionUDot
  NonlinearSystemBase & _nl;

  /// Name of the integrator's lumped mass diagonal vector
  const std::string _mass_vector_name;

  /// Computed lumped kinetic energy
  Real _ke;
};
