//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "AuxKernel.h"

/**
 * Sets the old solution value as the bound for the PETSc VI solvers on a
 * node-associated nonlinear variable. Companion to NodalConstantBounds: the dummy aux
 * variable does NOT need to match the finite element type of the bounded variable,
 * enabling families like BERNSTEIN on TET10 that the AuxKernel system rejects at
 * second order. Bounding Bernstein DoF coefficients by their old values enforces
 * pointwise field irreversibility, since Bernstein shape functions are non-negative.
 */
class NodalVariableOldValueBounds : public AuxKernel
{
public:
  static InputParameters validParams();

  NodalVariableOldValueBounds(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  enum BoundType
  {
    UPPER = 0,
    LOWER = 1
  };

  /// The type of bound (upper or lower)
  BoundType _type;

  /// The PETSc bound vector this object writes into
  NumericVector<Number> & _bounded_vector;

  /// The nonlinear variable being bounded
  MooseVariableFieldBase & _bounded_var;
};
