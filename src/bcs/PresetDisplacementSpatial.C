#include "PresetDisplacementSpatial.h"
#include "Function.h"

registerMooseObject("raccoonApp", PresetDisplacementSpatial);

InputParameters
PresetDisplacementSpatial::validParams()
{
  InputParameters params = DirichletBCBase::validParams();
  params.addClassDescription(
      "Newmark-consistent prescribed-displacement BC.  Evaluates the driving "
      "function's time derivative at the current node's position (or at "
      "coupled AuxVariable values via coupled_x / coupled_y / coupled_z, "
      "useful for restart simulations whose mesh is the dump-deformed config "
      "but whose undeformed coordinates are recovered as AuxVariables) so "
      "spatially-varying loads work correctly.  The stock PresetDisplacement "
      "evaluates the function at (0,0,0).");

  params.addParam<Real>("scale_factor", 1, "Scale factor if function is given.");
  params.addParam<FunctionName>("function", "1", "Function describing the displacement.");
  params.addRequiredCoupledVar("velocity", "The velocity variable.");
  params.addRequiredCoupledVar("acceleration", "The acceleration variable.");
  params.addRequiredParam<Real>("beta", "beta parameter for Newmark time integration.");

  params.addCoupledVar("coupled_x",
                        "Optional AuxVariable carrying the x-coordinate at which the "
                        "function should be evaluated.  Defaults to the node's own x.");
  params.addCoupledVar("coupled_y",
                        "Optional AuxVariable carrying the y-coordinate at which the "
                        "function should be evaluated.  Defaults to the node's own y.");
  params.addCoupledVar("coupled_z",
                        "Optional AuxVariable carrying the z-coordinate at which the "
                        "function should be evaluated.  Defaults to the node's own z.");

  // Forcefully preset the BC
  params.set<bool>("preset") = true;
  params.suppressParameter<bool>("preset");

  return params;
}

PresetDisplacementSpatial::PresetDisplacementSpatial(const InputParameters & parameters)
  : DirichletBCBase(parameters),
    _u_old(valueOld()),
    _scale_factor(parameters.get<Real>("scale_factor")),
    _function(getFunction("function")),
    _vel_old(coupledValueOld("velocity")),
    _accel_old(coupledValueOld("acceleration")),
    _beta(getParam<Real>("beta")),
    _has_coupled_x(isCoupled("coupled_x")),
    _has_coupled_y(isCoupled("coupled_y")),
    _has_coupled_z(isCoupled("coupled_z")),
    _coupled_x(_has_coupled_x ? &coupledValue("coupled_x") : nullptr),
    _coupled_y(_has_coupled_y ? &coupledValue("coupled_y") : nullptr),
    _coupled_z(_has_coupled_z ? &coupledValue("coupled_z") : nullptr)
{
}

Point
PresetDisplacementSpatial::spatialPoint() const
{
  const Real x = _has_coupled_x ? (*_coupled_x)[_qp] : (*_current_node)(0);
  const Real y = _has_coupled_y ? (*_coupled_y)[_qp] : (*_current_node)(1);
  const Real z = _has_coupled_z ? (*_coupled_z)[_qp] : (*_current_node)(2);
  return Point(x, y, z);
}

Real
PresetDisplacementSpatial::computeQpValue()
{
  // Evaluate the function's time derivative at the BC node's (potentially
  // overridden by AuxVariable) spatial position so spatially-varying loads
  // are honored.
  const Point p = spatialPoint();
  const Real vel = _function.timeDerivative(_t, p);
  const Real vel_old = _function.timeDerivative(_t - _dt, p);
  const Real accel = (vel - vel_old) / _dt;

  return _u_old[_qp] + _dt * _vel_old[_qp] +
         ((0.5 - _beta) * _accel_old[_qp] + _beta * accel) * _dt * _dt;
}
