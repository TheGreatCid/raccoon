
#pragma once

#include "AuxKernel.h"

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
  const VariableValue & _u_old_store;
  const VariableValue & _accel_old_store;
};
