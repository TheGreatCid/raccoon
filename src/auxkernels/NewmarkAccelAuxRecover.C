//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseTypes.h"
#include "NewmarkAccelAuxRecover.h"

registerMooseObject("raccoonApp", NewmarkAccelAuxRecover);

InputParameters
NewmarkAccelAuxRecover::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Computes the current acceleration using the Newmark method.");
  params.addRequiredCoupledVar("displacement", "displacement variable");
  params.addRequiredCoupledVar("velocity", "velocity variable");
  params.addRequiredParam<Real>("beta", "beta parameter for Newmark method");
  params.addRequiredParam<UserObjectName>("solution",
                                          "The SolutionUserObject to extract data from.");
  params.addRequiredParam<VariableName>("accel_name",
                                        "name of accel variable to get from solution object");

  return params;
}

NewmarkAccelAuxRecover::NewmarkAccelAuxRecover(const InputParameters & parameters)
  : AuxKernel(parameters),
    _disp_old(coupledValueOld("displacement")),
    _disp(coupledValue("displacement")),
    _vel_old(coupledValueOld("velocity")),
    _u_old(uOld()),
    _beta(getParam<Real>("beta")),
    _solution_object_ptr(NULL),
    _accel_name(getParam<VariableName>("accel_name"))
{
  _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");
}

Real
NewmarkAccelAuxRecover::computeValue()
{
  if (!isNodal())
    mooseError("must run on a nodal variable");
  Real accel_old;

  // Getting old acceleration from solution file
  if (_t_step == 0)
  {
    Point curr_Point = _q_point[_qp];
    accel_old = _solution_object_ptr->pointValue(1, curr_Point, _accel_name, nullptr);
  }
  else
    accel_old = _u_old[_qp];

  if (_t_step == 0)
    return accel_old;

  // Calculates acceeleration using Newmark time integration method
  return 1.0 / _beta *
         ((_disp[_qp] - _disp_old[_qp]) / (_dt * _dt) - _vel_old[_qp] / _dt -
          accel_old * (0.5 - _beta));
}
