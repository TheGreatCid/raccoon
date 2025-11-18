//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPenaltyConstraint.h"
#include "MooseUtils.h"
#include "RaccoonUtils.h"
#include "metaphysicl/raw_type.h"
registerMooseObject("raccoonApp", ADPenaltyConstraint);

InputParameters
ADPenaltyConstraint::validParams()
{
  InputParameters params = ADKernel::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("Enforces inequality constraints using the penalty method. "
                             "Can enforce lower bounds (irreversibility) or upper bounds on variables.");

  params.addParam<Real>(
      "penalty_param",
      1e4,
      "The penalty parameter. Larger values enforce constraints more strictly but hurt "
      "conditioning. Typical range: 1e2-1e6. Use adaptive_penalty=true for automatic tuning.");
  params.addParam<Real>(
      "epsilon",
      0.01,
      "Smoothing parameter for smooth approximations. Controls transition width. "
      "Larger values (0.01-0.1) improve convergence but relax constraint enforcement. "
      "Should scale with typical variable magnitudes.");
  params.addParam<bool>(
      "smooth", true, "Use smooth approximation of Macaulay bracket. Strongly recommended for "
                      "Newton convergence.");
  params.addParam<bool>("conditional", false, "Use conditional bounds?");
  params.addParam<bool>("upper", false, "Enforce upper bound (true) or lower bound (false)");
  params.addParam<Real>("upper_val", 1.0, "Upper bound value");

  MooseEnum smoothing_type("SQRT EXPONENTIAL POLYNOMIAL", "EXPONENTIAL");
  params.addParam<MooseEnum>("smoothing_type",
                             smoothing_type,
                             "Type of smooth approximation: "
                             "SQRT (current), EXPONENTIAL (better convergence), "
                             "POLYNOMIAL (C-infinity smooth)");

  params.addParam<bool>("adaptive_penalty",
                        false,
                        "Use adaptive penalty that increases if constraint is violated. "
                        "Improves convergence by starting with smaller penalty.");
  params.addParam<Real>("initial_penalty",
                        1e2,
                        "Initial penalty for adaptive scheme (adaptive_penalty=true)");
  params.addParam<Real>("penalty_growth",
                        2.0,
                        "Growth factor for adaptive penalty per Newton iteration where constraint "
                        "is violated significantly");
  params.addParam<Real>("violation_tol",
                        1e-3,
                        "Tolerance for constraint violation in adaptive scheme. "
                        "If violation > tol, increase penalty.");

  params.addParam<Real>("transition_width",
                        0.05,
                        "Width for smooth transitions in conditional mode (around d=0.95). "
                        "Larger values (0.05-0.1) improve convergence.");

  params.addParamNamesToGroup(
      "smoothing_type epsilon transition_width", "Smoothing and Regularization");
  params.addParamNamesToGroup(
      "adaptive_penalty initial_penalty penalty_growth violation_tol", "Adaptive Penalty");

  return params;
}

ADPenaltyConstraint::ADPenaltyConstraint(const InputParameters & parameters)
  : ADKernel(parameters),
    BaseNameInterface(parameters),
    _penalty(getParam<Real>("penalty_param")),
    _epsilon(getParam<Real>("epsilon")),
    _u_old(valueOld()),
    _smooth(getParam<bool>("smooth")),
    _conditional(getParam<bool>("conditional")),
    _upper(getParam<bool>("upper")),
    _upper_val(getParam<Real>("upper_val")),
    _smoothing_type(getParam<MooseEnum>("smoothing_type")),
    _adaptive_penalty(getParam<bool>("adaptive_penalty")),
    _initial_penalty(getParam<Real>("initial_penalty")),
    _penalty_growth(getParam<Real>("penalty_growth")),
    _violation_tol(getParam<Real>("violation_tol")),
    _transition_width(getParam<Real>("transition_width"))
{
  if (_adaptive_penalty && _penalty != 1e4)
    mooseWarning("ADPenaltyConstraint: adaptive_penalty is enabled but penalty_param is set. "
                 "The penalty_param will be used as initial value, but initial_penalty is "
                 "recommended for adaptive mode.");

  if (!_smooth && (_smoothing_type != "SQRT"))
    mooseWarning("ADPenaltyConstraint: smoothing_type is set but smooth=false. "
                 "Set smooth=true to use the specified smoothing method.");

  if (_epsilon < 1e-6)
    mooseWarning("ADPenaltyConstraint: epsilon=",
                 _epsilon,
                 " is very small. This creates sharp transitions and may hurt convergence. "
                 "Consider using epsilon in range [0.01, 0.1] for better Newton convergence.");

  if (_penalty > 1e6)
    mooseWarning("ADPenaltyConstraint: penalty_param=",
                 _penalty,
                 " is very large. This creates ill-conditioned systems and may cause divergence. "
                 "Consider using adaptive_penalty=true or smaller penalty with larger epsilon.");
}

ADReal
ADPenaltyConstraint::smoothHeaviside(const ADReal x, const Real width) const
{
  // Smooth approximation: H(x) ≈ 0.5 * (1 + tanh(x/width))
  // Continuous and differentiable, transitions over ~4*width
  return 0.5 * (1.0 + std::tanh(x / width));
}

ADReal
ADPenaltyConstraint::smoothMacaulay(const ADReal x, const Real eps) const
{
  if (_smoothing_type == "SQRT")
  {
    // Original: f(x) = 0.5 * (sqrt(x² + ε²) + x)
    return 0.5 * (std::sqrt(x * x + eps * eps) + x);
  }
  else if (_smoothing_type == "EXPONENTIAL")
  {
    // Exponential smoothing: f(x) = x + ε*ln(1 + exp(-x/ε))
    // Much better conditioned than SQRT, smoother derivatives
    // For large |x/ε|, approaches max(0,x) exponentially fast
    const ADReal z = -x / eps;
    if (MetaPhysicL::raw_value(z) > 20.0) // exp(-z) ≈ 0, avoid overflow
      return x > 0 ? x : eps * z; // Returns x for x>0, nearly 0 for x<0
    else if (MetaPhysicL::raw_value(z) < -20.0) // exp(-z) huge, but log cancels
      return eps * std::log(std::exp(z));
    else
      return x + eps * std::log(1.0 + std::exp(z));
  }
  else // POLYNOMIAL
  {
    // Polynomial smoothing (C-infinity): quintic interpolation
    // f(x) = 0 for x < -ε, x for x > ε, smooth cubic in between
    if (x < -eps)
      return 0.0;
    else if (x > eps)
      return x;
    else
    {
      // Quintic: f(t) = t + ε*(6t⁵ - 15t⁴ + 10t³)/16 where t = (x+ε)/(2ε) ∈ [0,1]
      const ADReal t = (x + eps) / (2.0 * eps);
      const ADReal t3 = t * t * t;
      const ADReal t4 = t3 * t;
      const ADReal t5 = t4 * t;
      return t * (2.0 * eps) - eps + eps * (6.0 * t5 - 15.0 * t4 + 10.0 * t3);
    }
  }
}

ADReal
ADPenaltyConstraint::computeQpResidual()
{
  ADReal function = 0;
  ADReal delta_d;

  if (!_upper)
  {
    // Lower bound constraint: enforce u >= u_old (or u >= 0)
    if (!_conditional)
    {
      // Smooth version of: delta_d = (u_old > 0) ? (u_old - u) : (0 - u)
      const ADReal u_old_positive = smoothHeaviside(_u_old[_qp], _epsilon);
      delta_d = u_old_positive * (_u_old[_qp] - _u[_qp]) + (1.0 - u_old_positive) * (0 - _u[_qp]);
    }
    else
    {
      // Smooth version of: if (u_old > 0.95) use u_old-u, else use 0-u
      // Smooth transition centered at 0.95 with width _transition_width
      const ADReal weight = smoothHeaviside(_u_old[_qp] - 0.95, _transition_width);
      delta_d = weight * (_u_old[_qp] - _u[_qp]) + (1.0 - weight) * (0 - _u[_qp]);
    }
  }
  else
  {
    // Upper bound constraint: enforce u <= upper_val
    delta_d = _u[_qp] - _upper_val;
  }

  // Apply smooth Macaulay bracket to constraint violation
  if (_smooth)
    function = smoothMacaulay(delta_d, _epsilon);
  else
    function = RaccoonUtils::Macaulay(delta_d);

  // Adaptive penalty: increase penalty if constraint is significantly violated
  Real penalty = _penalty;
  if (_adaptive_penalty)
  {
    const Real violation = MetaPhysicL::raw_value(function);
    if (violation > _violation_tol)
    {
      // Increase penalty logarithmically with violation magnitude
      const Real violation_factor = std::log10(std::max(violation / _violation_tol, 1.0)) + 1.0;
      penalty = _initial_penalty * std::pow(_penalty_growth, violation_factor);
      penalty = std::min(penalty, 1e8); // Cap to prevent overflow
    }
    else
    {
      penalty = _initial_penalty;
    }
  }

  const Real sign = _upper ? +1.0 : -1.0; // upper‑bound uses +, lower‑bound uses –
  return sign * penalty * _test[_i][_qp] * function;
}
