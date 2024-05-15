//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NewmarkVelAuxRecover.h"

registerMooseObject("raccoonApp", NewmarkVelAuxRecover);

InputParameters
NewmarkVelAuxRecover::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Calculates the current velocity using Newmark method.");
  params.addRequiredCoupledVar("acceleration", "acceleration variable");
  params.addRequiredParam<Real>("gamma", "gamma parameter for Newmark method");
  params.addRequiredParam<UserObjectName>("solution",
                                          "The SolutionUserObject to extract data from.");
  params.addRequiredParam<VariableName>("vel_name",
                                        "name of vel variable to get from solution object");
  return params;
}

NewmarkVelAuxRecover::NewmarkVelAuxRecover(const InputParameters & parameters)
  : AuxKernel(parameters),
    _accel_old(coupledValueOld("acceleration")),
    _accel(coupledValue("acceleration")),
    _u_old(uOld()),
    _gamma(getParam<Real>("gamma")),
    _solution_object_ptr(NULL),
    _vel_name(getParam<VariableName>("vel_name"))
{
  _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");
}

Real
NewmarkVelAuxRecover::computeValue()
{
  Real vel_old;
  if (_t_step == 0)
  {
    Point curr_Point = _q_point[_qp];
    vel_old = _solution_object_ptr->pointValue(1, curr_Point, _vel_name, nullptr);
  }
  else
  {
    vel_old = _u_old[_qp];
  }
  if (!isNodal())
    mooseError("must run on a nodal variable");
  // Calculates Velocity using Newmark time integration scheme
  if (_t_step == 0)
  {
    return vel_old;
  }
  else
  {
    return vel_old + (_dt * (1 - _gamma)) * _accel_old[_qp] + _gamma * _dt * _accel[_qp];
  }
}
