#include "ADInertialForceRecover.h"
#include "SubProblem.h"
#include "TimeIntegrator.h"
#include "NonlinearSystemBase.h"

registerMooseObject("raccoonApp", ADInertialForceRecover);

InputParameters
ADInertialForceRecover::validParams()
{
  InputParameters params = ADTimeKernel::validParams();
  params.addClassDescription("Calculates the residual for the inertial force "
                             "($M \\cdot acceleration$) and the contribution of mass"
                             " dependent Rayleigh damping and HHT time "
                             " integration scheme ($\\eta \\cdot M \\cdot"
                             " ((1+\\alpha)velq2-\\alpha \\cdot vel-old) $)");
  params.set<bool>("use_displaced_mesh") = true;
  params.addCoupledVar("velocity", "velocity variable");
  params.addCoupledVar("acceleration", "acceleration variable");
  params.addParam<Real>("beta", "beta parameter for Newmark Time integration");
  params.addParam<Real>("gamma", "gamma parameter for Newmark Time integration");
  params.addParam<MaterialPropertyName>("eta",
                                        0.0,
                                        "Name of material property or a constant real "
                                        "number defining the eta parameter for the "
                                        "Rayleigh damping.");
  params.addParam<MaterialPropertyName>(
      "density_scaling",
      0.0,
      "Name of material property to add mass scaling in explicit simulations.");
  params.addParam<Real>("alpha",
                        0,
                        "alpha parameter for mass dependent numerical damping induced "
                        "by HHT time integration scheme");
  params.addParam<MaterialPropertyName>(
      "density", "density", "Name of Material Property that provides the density");
  params.addRequiredParam<VariableName>("inert_name",
                                        "name of inertial variable to get from solution object");
  return params;
}

ADInertialForceRecover::ADInertialForceRecover(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _density(this->getGenericMaterialProperty<Real, true>("density")),
    _has_beta(this->isParamValid("beta")),
    _has_gamma(this->isParamValid("gamma")),
    _beta(_has_beta ? this->getParam<Real>("beta") : 0.1),
    _gamma(_has_gamma ? this->getParam<Real>("gamma") : 0.1),
    _has_velocity(this->isParamValid("velocity")),
    _has_acceleration(this->isParamValid("acceleration")),
    _eta(this->getGenericMaterialProperty<Real, true>("eta")),
    _density_scaling(this->getMaterialProperty<Real>("density_scaling")),
    _alpha(this->getParam<Real>("alpha")),
    _time_integrator(*_sys.getTimeIntegrator()),
    _solution_object_ptr(NULL),
    _inert_name(getParam<VariableName>("inert_name"))

{
  if (_has_beta && _has_gamma && _has_velocity && _has_acceleration)
  {
    _vel_old = &this->coupledValueOld("velocity");
    _accel_old = &this->coupledValueOld("acceleration");
    _u_old = &this->valueOld();
  }
  else if (!_has_beta && !_has_gamma && !_has_velocity && !_has_acceleration)
  {
    _u_dot_old = &(_var.uDotOld());
    _du_dot_du = &(_var.duDotDu());
    _du_dotdot_du = &(_var.duDotDotDu());

    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotFactorTag());
    this->addFEVariableCoupleableVectorTag(_time_integrator.uDotDotFactorTag());

    _u_dot_factor = &_var.vectorTagValue(_time_integrator.uDotFactorTag());
    _u_dotdot_factor = &_var.vectorTagValue(_time_integrator.uDotDotFactorTag());

    if (_time_integrator.isLumped())
    {
      _u_dot_factor_dof = &_var.vectorTagDofValue(_time_integrator.uDotFactorTag());
      _u_dotdot_factor_dof = &_var.vectorTagDofValue(_time_integrator.uDotDotFactorTag());
    }
  }
  else
    mooseError("InertialForce: Either all or none of `beta`, `gamma`, `velocity`and `acceleration` "
               "should be provided as input.");

  // Check if HHT and explicit are being used simultaneously
  if (_alpha != 0 && _time_integrator.isExplicit())
    mooseError("InertialForce: HHT time integration parameter can only be used with Newmark-Beta "
               "time integration.");

  // Check if beta and explicit are being used simultaneously
  if (_has_beta && _time_integrator.isExplicit())
    mooseError("InertialForce: Newmark-beta integration parameter, beta, cannot be provided along "
               "with an explicit time "
               "integrator.");

  if (!(_time_integrator.isLumped() && _time_integrator.isExplicit()))
    if (parameters.isParamSetByUser("density_scaling"))
      this->paramError(
          "density_scaling",
          "Density (mass) scaling can only be used in lumped mass, explicit simulations");
}

GenericReal<true>
ADInertialForceRecover::computeQpResidual()
{

  if (_dt == 0)
    return 0;
  else if (_has_beta)
  {
    if (_t_step == 1)
    {
      Point curr_Point = _q_point[_qp];
      ADReal inert_old = _solution_object_ptr->pointValue(1, curr_Point, _inert_name, nullptr);

      auto accel = 1.0 / _beta *
                   (((_u[_qp] - inert_old) / (_dt * _dt)) - (*_vel_old)[_qp] / _dt -
                    (*_accel_old)[_qp] * (0.5 - _beta));
      auto vel =
          (*_vel_old)[_qp] + (_dt * (1.0 - _gamma)) * (*_accel_old)[_qp] + _gamma * _dt * accel;
      return _test[_i][_qp] * _density[_qp] *
             (accel + vel * _eta[_qp] * (1.0 + _alpha) - _alpha * _eta[_qp] * (*_vel_old)[_qp]);
    }
    else
    {
      auto accel = 1.0 / _beta *
                   (((_u[_qp] - (*_u_old)[_qp]) / (_dt * _dt)) - (*_vel_old)[_qp] / _dt -
                    (*_accel_old)[_qp] * (0.5 - _beta));
      auto vel =
          (*_vel_old)[_qp] + (_dt * (1.0 - _gamma)) * (*_accel_old)[_qp] + _gamma * _dt * accel;
      return _test[_i][_qp] * _density[_qp] *
             (accel + vel * _eta[_qp] * (1.0 + _alpha) - _alpha * _eta[_qp] * (*_vel_old)[_qp]);
    }
  }
  // Lumped mass option
  // Only lumping the masses here
  // will multiply by corresponding residual multiplier after lumping the matrix
  // Density scaling is a fictitious added density to increase the time step
  else if (_time_integrator.isLumped() && _time_integrator.isExplicit() && true)
    return _test[_i][_qp] * (_density[_qp] + _density_scaling[_qp]);

  // Consistent mass option
  // Same for explicit, implicit, and implicit with HHT
  else
    return _test[_i][_qp] * _density[_qp] *
           ((*_u_dotdot_factor)[_qp] + (*_u_dot_factor)[_qp] * _eta[_qp] * (1.0 + _alpha) -
            _alpha * _eta[_qp] * (*_u_dot_old)[_qp]);
}

void
ADInertialForceRecover::computeResidualAdditional()
{
  mooseError("Internal error");
}

Real
ADInertialForceRecover::computeQpJacobian()
{
  mooseError("Internal error");
}
