//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NewmarkVelAuxRecover.h"

registerMooseObject("RACCOON", NewmarkVelAuxRecover);

InputParameters
NewmarkVelAuxRecover::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the current velocity using Newmark method.");
  params.addRequiredCoupledVar("acceleration", "acceleration variable");
  params.addRequiredParam<Real>("gamma", "gamma parameter for Newmark method");
  params.addRequiredCoupledVar("u_old_store", "u_old_store");
  params.addRequiredCoupledVar("accel_old_store", "accel_old_store");

  return params;
}

NewmarkVelAuxRecover::NewmarkVelAuxRecover(const InputParameters & parameters)
  : AuxKernel(parameters),
    _accel_old(coupledValueOld("acceleration")),
    _accel(coupledValue("acceleration")),
    _u_old(uOld()),
    _gamma(getParam<Real>("gamma")),
    _u_old_store(coupledValue("u_old_store")),
    _accel_old_store(coupledValue("accel_old_store"))
{
}

Real
NewmarkVelAuxRecover::computeValue()
{
  Real vel_old;
  if (_t_step == 1)
  {
    vel_old = _u_old_store[_qp];
  }
  else
  {
    vel_old = _u_old[_qp];
  }
  if (!isNodal())
    mooseError("must run on a nodal variable");
  // Calculates Velocity using Newmark time integration scheme
  if (_t_step == 1)
  {
    return vel_old + (_dt * (1 - _gamma)) * _accel_old_store[_qp] + _gamma * _dt * _accel[_qp];
  }
  else
  {
    return vel_old + (_dt * (1 - _gamma)) * _accel_old[_qp] + _gamma * _dt * _accel[_qp];
  }
}
