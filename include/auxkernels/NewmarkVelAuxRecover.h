
#pragma once

#include "AuxKernel.h"
#include "SolutionUserObject.h"

class NewmarkVelAuxRecover : public AuxKernel
{
public:
  static InputParameters validParams();

  /**
   * Calcualtes velocity using Newmark time integration scheme
   */
  NewmarkVelAuxRecover(const InputParameters & parameters);

  virtual ~NewmarkVelAuxRecover() {}

protected:
  virtual Real computeValue();

  const VariableValue & _accel_old;
  const VariableValue & _accel;
  const VariableValue & _u_old;
  Real _gamma;
  const SolutionUserObject * _solution_object_ptr;
  VariableName _vel_name;
};
