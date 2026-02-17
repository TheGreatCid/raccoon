//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "JohnsonCookHardening.h"

registerMooseObject("raccoonApp", JohnsonCookHardening);

InputParameters
JohnsonCookHardening::validParams()
{
  InputParameters params = PlasticHardeningModel::validParams();
  params.addClassDescription("The Johnson-Cook plasticity model.");

  params.addRequiredParam<MaterialPropertyName>("sigma_0",
                                                "The reference yield stress $\\sigma_0$");
  params.addRequiredParam<Real>("n", "The exponent n in the JC Model");
  params.addRequiredParam<Real>("m", "The exponent m in the JC Model");
  params.addRequiredParam<Real>("reference_plastic_strain",
                                "The $\\epsilon_0$ parameter in the JC Model");
  params.addRequiredParam<Real>("reference_plastic_strain_rate",
                                "The ref plastic strain rate parameter in the JC model");
  params.addRequiredParam<Real>("T0", "Reference temperature of the material");
  params.addRequiredCoupledVar("T", "Temperature");
  params.addRangeCheckedParam<Real>(
      "taylor_quinney_factor",
      1,
      "taylor_quinney_factor<=1 & taylor_quinney_factor>=0",
      "The Taylor-Quinney factor. 1 (default) for purely dissipative, 0 for purely energetic.");
  params.addRequiredParam<MaterialPropertyName>("A", "'A' parameter for the JC model");
  params.addRequiredParam<Real>("B", "'B' parameter for the JC model");
  params.addRequiredParam<Real>("C", "'C' parameter for the JC model");
  params.addRequiredParam<Real>("Tm", "The melting temperature of the material");

  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "plastic_energy_density",
      "psip",
      "Name of the plastic energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "gp", "The degradation function");
  params.addParam<bool>("disable_dissipation",
                        false,
                        "Set to true to turn off Johnson-Cook plastic dissipation contributions.");
  return params;
}

JohnsonCookHardening::JohnsonCookHardening(const InputParameters & parameters)
  : PlasticHardeningModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _sigma_0(getADMaterialProperty<Real>(prependBaseName("sigma_0", true))),
    _n(getParam<Real>("n")),
    _m(getParam<Real>("m")),
    _ep0(getParam<Real>("reference_plastic_strain")),
    _epdot0(getParam<Real>("reference_plastic_strain_rate")),
    _T0(getParam<Real>("T0")),
    _T(coupledValue("T")),
    _T_old(coupledValueOld("T")),
    _tqf(getParam<Real>("taylor_quinney_factor")),
    _A(getADMaterialProperty<Real>("A")),
    _B(getParam<Real>("B")),
    _C(getParam<Real>("C")),
    _Tm(getParam<Real>("Tm")),
    _d_name(getVar("phase_field", 0)->name()),

    // The strain energy density and its derivatives
    _psip_name(prependBaseName("plastic_energy_density", true)),
    _psip(declareADProperty<Real>(_psip_name)),
    _psip_active(declareADProperty<Real>(_psip_name + "_active")),
    _dpsip_dd(declareADProperty<Real>(derivativePropertyName(_psip_name, {_d_name}))),

    // The degradation function and its derivatives
    _gp_name(prependBaseName("degradation_function", true)),
    _gp(getADMaterialProperty<Real>(_gp_name)),
    _dgp_dd(getADMaterialProperty<Real>(derivativePropertyName(_gp_name, {_d_name}))),
    _disable_dissipation(getParam<bool>("disable_dissipation")),
    _T_local(),
    _use_local_T()
{
  _T_local.resize(_fe_problem.getMaxQps());
  _use_local_T.resize(_fe_problem.getMaxQps(), 0);
}

ADReal // Temperature degradation
JohnsonCookHardening::temperatureDependence()
{
  using std::pow;
  return 1 - pow((getQpT() - _T0) / (_Tm - _T0), _m);
}

ADReal
JohnsonCookHardening::initialGuess(const ADReal & effective_trial_stress)
{
  using std::max;
  using std::pow;
  ADReal trial_over_stress =
      effective_trial_stress / _sigma_0[_qp] / temperatureDependence() - _A[_qp];
  if (trial_over_stress < 0)
    trial_over_stress = 0;
  return max(_ep0 * pow(trial_over_stress / _B, 1 / _n), libMesh::TOLERANCE * libMesh::TOLERANCE);
}

ADReal
JohnsonCookHardening::plasticEnergy(const ADReal & ep, const unsigned int derivative)
{
  if (derivative == 0)
  {
    using std::pow;
    _psip_active[_qp] = (1 - _tqf) * _sigma_0[_qp] *
                        (_A[_qp] * ep + _B * _ep0 * pow(ep / _ep0, _n + 1) / (_n + 1)) *
                        temperatureDependence();
    _psip[_qp] = _gp[_qp] * _psip_active[_qp];
    _dpsip_dd[_qp] = _dgp_dd[_qp] * _psip_active[_qp];
    return _psip[_qp];
  }

  if (derivative == 1)
  {
    using std::pow;
    return _gp[_qp] * (1 - _tqf) * _sigma_0[_qp] * (_A[_qp] + _B * pow(ep / _ep0, _n)) *
           temperatureDependence();
  }
  if (derivative == 2)
  {
    using std::pow;
    return _gp[_qp] * (1 - _tqf) * _sigma_0[_qp] * _B * pow(ep / _ep0, _n - 1) * _n / _ep0 *
           temperatureDependence();
  }
  mooseError(name(), "internal error: unsupported derivative order.");
}

ADReal
JohnsonCookHardening::plasticDissipation(const ADReal & delta_ep,
                                         const ADReal & ep,
                                         const unsigned int derivative)
{
  if (_disable_dissipation)
    return ADReal(0);

  // For all cases, we are splitting between rate dependent and non rate dependent portions to avoid
  // /0 errors

  ADReal result = 0;

  if (derivative == 0)
  {
    using std::log;
    using std::pow;
    result += (_A[_qp] + _B * pow(ep / _ep0, _n)) * _tqf * delta_ep;
    if (_t_step > 0 && delta_ep > libMesh::TOLERANCE * libMesh::TOLERANCE)
      result += (_A[_qp] + _B * pow(ep / _ep0, _n)) * (_C * log(delta_ep / _dt / _epdot0) - _C) *
                delta_ep;
  }

  if (derivative == 1)
  {
    using std::log;
    using std::pow;
    result += (_A[_qp] + _B * pow(ep / _ep0, _n)) * _tqf;
    if (_t_step > 0 && delta_ep > libMesh::TOLERANCE * libMesh::TOLERANCE)
      result += (_A[_qp] + _B * pow(ep / _ep0, _n)) * (_C * log(delta_ep / _dt / _epdot0));
  }

  if (derivative == 2)
  {
    using std::log;
    using std::pow;
    result += _B * pow(ep / _ep0, _n - 1) * _n / _ep0 * _tqf;
    if (_t_step > 0 && delta_ep > libMesh::TOLERANCE * libMesh::TOLERANCE)
      result += (_A[_qp] + _B * pow(ep / _ep0, _n)) * _C / delta_ep +
                _B * pow(ep / _ep0, _n - 1) * _n / _ep0 * _C * log(delta_ep / _dt / _epdot0);
  }

  return _gp[_qp] * result * _sigma_0[_qp] * temperatureDependence();

  mooseError(name(), "internal error: unsupported derivative order.");
}

Real
JohnsonCookHardening::temperatureDependenceLogDerivative(Real T)
{
  // xi = d(TD)/dT / TD = -m * (1-TD) / (TD * (T-T0))
  // Guards against T <= T0 (cold limit) and T >= Tm (melt).
  if (T <= _T0 + 1e-14 || T >= _Tm)
    return 0.0;
  using std::pow;
  const Real u = (T - _T0) / (_Tm - _T0);
  const Real TD = 1.0 - pow(u, _m);
  if (TD <= 1e-14)
    return 0.0;
  return -_m * (1.0 - TD) / (TD * (T - _T0));
}

Real
JohnsonCookHardening::thermalConjugateTemperatureDerivative(Real ep)
{
  // dTC/dT = TC * (m*T - T0) / (T * (T - T0))
  if (_disable_dissipation)
    return 0.0;
  const Real T = getQpT();
  if (T <= _T0 + 1e-14 || T >= _Tm)
    return 0.0;
  using std::pow;
  const Real u = (T - _T0) / (_Tm - _T0);
  const Real one_minus_TD = pow(u, _m);
  const Real A_ep = MetaPhysicL::raw_value(_A[_qp]) + _B * pow(ep / _ep0, _n);
  // TC from thermalConjugate, computed in Real arithmetic:
  // TC = gp * T * (1-tqf) * sigma_0 * A_ep * m * (1-TD) / (T0-T)
  const Real TC = MetaPhysicL::raw_value(_gp[_qp]) * T * (1.0 - _tqf) *
                  MetaPhysicL::raw_value(_sigma_0[_qp]) * A_ep * _m * one_minus_TD / (_T0 - T);
  return TC * (_m * T - _T0) / (T * (T - _T0));
}

Real
JohnsonCookHardening::dissipationFlowStressRateJacobian(Real dep, Real ep)
{
  // dep * partial(W_d_flowstress)/partial(dep)
  // Only the rate-dependent (log) term contributes: gp * sigma_0 * A_ep * C * TD
  if (_disable_dissipation || _t_step == 0 || dep <= libMesh::TOLERANCE * libMesh::TOLERANCE)
    return 0.0;
  using std::pow;
  const Real T = getQpT();
  if (T <= _T0 + 1e-14 || T >= _Tm)
    return 0.0;
  const Real TD = 1.0 - pow((T - _T0) / (_Tm - _T0), _m);
  const Real A_ep = MetaPhysicL::raw_value(_A[_qp]) + _B * pow(ep / _ep0, _n);
  return MetaPhysicL::raw_value(_gp[_qp]) * MetaPhysicL::raw_value(_sigma_0[_qp]) * A_ep * _C * TD;
}

ADReal // Thermal conjugate term
JohnsonCookHardening::thermalConjugate(const ADReal & ep)
{
  using std::pow;

  if (_disable_dissipation)
    return ADReal(0);
  const Real T = getQpT();
  return _gp[_qp] * T * (1 - _tqf) * _sigma_0[_qp] * (_A[_qp] + _B * pow(ep / _ep0, _n)) *
         (_m * (pow((_T0 - T) / (_T0 - _Tm), _m))) / (_T0 - T);
}

void
JohnsonCookHardening::setLocalTemperature(Real T)
{
  if (_T_local.size() == 0)
  {
    _T_local.resize(_fe_problem.getMaxQps());
    _use_local_T.resize(_fe_problem.getMaxQps(), 0);
  }
  _T_local[_qp] = T;
  _use_local_T[_qp] = 1;
}

void
JohnsonCookHardening::clearLocalTemperature()
{
  // assume sized already
  if (_use_local_T.size() > 0)
    _use_local_T[_qp] = 0;
}

Real
JohnsonCookHardening::getQpTemperatureOld() const
{
  // return previous-step temperature at this qp (used in some places)
  return _T_old[_qp];
}

Real
JohnsonCookHardening::getQpT() const
{
  if (_use_local_T.size() > 0 && _use_local_T[_qp])
    return _T_local[_qp];
  return _T[_qp]; // current iterative temperature DOF
}
