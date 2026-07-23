//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "AuxKernel.h"

/**
 * Sets a constant bound for the PETSc VI solvers on a node-associated nonlinear
 * variable. Unlike the framework's ConstantBounds, the dummy aux variable does NOT
 * need to match the finite element type of the bounded variable -- it only needs to
 * be a nodal (Lagrange) variable of matching order so that every node carrying a DoF
 * of the bounded variable is visited. This enables bounding variables whose families
 * are rejected by the AuxKernel system at second order, e.g. BERNSTEIN on TET10,
 * where each node carries exactly one DoF.
 */
class NodalConstantBounds : public AuxKernel
{
public:
  static InputParameters validParams();

  NodalConstantBounds(const InputParameters & parameters);

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

  /// The constant bound value
  const Real _bound_value;
};
