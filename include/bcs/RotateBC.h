//* This file is part of the RACCOON framework

#pragma once

#include "ADDirichletBCBase.h"
/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the values of a nodal variable at nodes
 */
class RotateBC : public ADDirichletBCBase
{
public:
  static InputParameters validParams();

  RotateBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpValue() override;
  const VariableValue & _theta;

  /// the displacement value
  MooseEnum _xy;
};
