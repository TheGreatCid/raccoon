//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "AuxKernel.h"

/**
 * Copies the time derivative (`coupledDot`) of a coupled variable into an
 * AuxVariable.  Works for both elemental and nodal AuxVariables (unlike
 * MOOSE's `TimeDerivativeAux`, which is elemental-only).
 *
 * Use case: with the ExplicitMixedOrder / central-difference time integrator,
 * the variable's time derivative (`_var.uDot()`) is the integrator's internal
 * velocity vector.  Mirroring that to a nodal AuxVariable lets it be dumped
 * to exodus for downstream recovery / postprocessing.
 *
 * Set `order_derivative = SECOND` to instead copy the second time derivative
 * (`coupledDotDot`, the integrator's acceleration for CD).
 */
class CoupledTimeDerivativeAux : public AuxKernel
{
public:
  static InputParameters validParams();
  CoupledTimeDerivativeAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// Time derivative (or second time derivative if `order_derivative = SECOND`)
  /// of the coupled variable.  Indexed by [_qp]; for nodal aux _qp is always 0.
  const VariableValue & _coupled_dot;
};
