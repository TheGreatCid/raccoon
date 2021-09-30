//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "PowerLawHardening.h"
#include <vector>
registerMooseObject("raccoonApp", PowerLawHardening);

InputParameters
PowerLawHardening::validParams()
{
  InputParameters params = PlasticHardeningModel::validParams();
  params.addClassDescription("Plastic hardening following a power law.");

  params.addRequiredParam<MaterialPropertyName>("exponent",
                                                "The exponent n in the power law hardening $n$");
  params.addRequiredParam<MaterialPropertyName>(
      "reference_plastic_strain", "The $\\epsilon_0$ parameter in the power law hardening");

  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "plastic_energy_density",
      "psip",
      "Name of the plastic energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "gp", "The degradation function");
  params.addRequiredParam<MaterialPropertyName>("ref_yield_stress",
                                                "The reference yield stress $\\sigma_0$");
  params.addRequiredParam<MaterialPropertyName>("T0", "T0");
  params.addParam<MooseEnum>(
      "sigy_func", MooseEnum("PIECE EXP", "PIECE"), "The function for degrading yield stress");

  return params;
}

PowerLawHardening::PowerLawHardening(const InputParameters & parameters)
  : PlasticHardeningModel(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _n(getADMaterialProperty<Real>(prependBaseName("exponent", true))),
    _ep0(getADMaterialProperty<Real>(prependBaseName("reference_plastic_strain", true))),

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
    _T(getMaterialProperty<Real>(prependBaseName("Temp"))), // Remove
    _sigma_0(getADMaterialProperty<Real>(prependBaseName("ref_yield_stress"))),
    _sigma_y(declareADProperty<Real>(prependBaseName("yield_stress"))),
    _T0(getADMaterialProperty<Real>(prependBaseName("T0"))),
    _sigy_func(getParam<MooseEnum>("sigy_func").getEnum<Sigy_func>())
{
}
ADReal
PowerLawHardening::plasticEnergy(const ADReal & ep, const unsigned int derivative)
{
  if (_sigy_func == Sigy_func::exp)
  {
    _sigma_y[_qp] = _sigma_0[_qp] * (std::exp((_T0[_qp] - _T[_qp]) / 500));
  }
  else if (_sigy_func == Sigy_func::piece)
  {
    _sigma_y[_qp] = _sigma_0[_qp] * piecewise(); //(std::exp((_T0[_qp] - _T[_qp]) / 500));
  }
  else
  {
    paramError("sigy_func", "Wrong function input");
  }
  if (derivative == 0)
  {
    _psip_active[_qp] = _n[_qp] * _sigma_y[_qp] * _ep0[_qp] / (_n[_qp] + 1) *
                        (std::pow(1 + ep / _ep0[_qp], 1 / _n[_qp] + 1) - 1);
    _psip[_qp] = _gp[_qp] * _psip_active[_qp];
    _dpsip_dd[_qp] = _dgp_dd[_qp] * _psip_active[_qp];
    return _psip[_qp];
  }

  if (derivative == 1)
    return _gp[_qp] * _sigma_y[_qp] * std::pow(1 + ep / _ep0[_qp], 1 / _n[_qp]);

  if (derivative == 2)
    return _gp[_qp] * _sigma_y[_qp] * std::pow(1 + ep / _ep0[_qp], 1 / _n[_qp] - 1) / _n[_qp] /
           _ep0[_qp];

  mooseError(name(), "internal error: unsupported derivative order.");
  return 0;
}

ADReal
PowerLawHardening::piecewise()
{

  Real Tl = _T[_qp];

  if (273.15 <= Tl && Tl < 332.9)
  {
    std::vector<Real> coeff = {
        -5.61270528573491e-08, 1.58495556827408e-05, -0.00158192359946669, 1};
    return val(Tl, coeff, 273.15);
  }
  else if (332.9 <= Tl && Tl < 394.3)
  {
    std::vector<Real> coeff = {
        -5.61270528573493e-08, 5.78746500533272e-06, -0.000288942448851826, 0.950089179000000};
    return val(Tl, coeff, 332.9);
  }
  else if (394.3 <= Tl && Tl < 512.6)
  {
    std::vector<Real> coeff = {
        2.22603215416878e-08, -4.54210463506049e-06, -0.000212544134567699, 0.941186000000000};
    return val(Tl, coeff, 394.3);
  }
  else if (512.6 <= Tl && Tl < 684.4)
  {
    std::vector<Real> coeff = {
        -1.13978549501960e-08, 3.36085514388411e-06, -0.000352334975716497, 0.889315233000000};
    return val(Tl, coeff, 512.6);
  }
  else if (684.4 <= Tl && Tl < 813.5)
  {
    std::vector<Real> coeff = {
        -2.44797253852972e-08, -2.51268543339399e-06, -0.000206642087807611, 0.870190744000000};
    return val(Tl, coeff, 684.4);
  }
  else if (813.5 <= Tl && Tl < 873.3)
  {
    std::vector<Real> coeff = {
        8.58904452147939e-08, -1.19956310052036e-05, -0.00208005056441911, 0.748906974000000};
    return val(Tl, coeff, 813.5);
  }
  else if (873.3 <= Tl && Tl < 933.8)
  {
    std::vector<Real> coeff = {
        -1.17165095137872e-07, 3.41789111374217e-06, -0.00259315840863068, 0.599942456000000};
    return val(Tl, coeff, 873.3);
  }
  else if (933.841425300000 <= Tl && Tl < 993.616700000000)
  {
    std::vector<Real> coeff = {
        2.34876403207678e-07, -1.78572109163676e-05, -0.00346713315251594, 0.429526021000000};
    return val(Tl, coeff, 933.841425300000);
  }
  else if (993.61670000000 <= Tl && Tl < 1103.3)
  {
    std::vector<Real> coeff = {
        -6.15071218938525e-08, 2.42621936504931e-05, -0.00308427355013482, 0.208637262000000};
    return val(Tl, coeff, 993.61670000000);
  }
  else if (1103.3 <= Tl && Tl < 1213.7)
  {
    std::vector<Real> coeff = {
        -6.15071218938525e-08, 4.03227602052467e-06, 1.77780129137437e-05, 0.0810667030000000};
    return val(Tl, coeff, 1103.3);
  }
  // std::cout << _T[_qp] << std::endl;
  // std::cout << Tl << std::endl;
  //  paramError("T0", "outside threshhold");

  return 1;
}

Real
PowerLawHardening::val(Real Tl, const std::vector<Real> & coeff, Real breaks)
{
  return coeff[0] * std::pow((Tl - breaks), (4 - 1)) + coeff[1] * std::pow((Tl - breaks), (4 - 2)) +
         coeff[2] * std::pow((Tl - breaks), (4 - 3)) + coeff[3];
}
