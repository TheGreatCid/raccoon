//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "PsicDeg.h"
registerMooseObject("raccoonApp", PsicDeg);

InputParameters
PsicDeg::validParams()
{
  InputParameters params = Material::validParams();

  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");
  params.addRequiredParam<Real>("A", "'A' parameter for the psic func");
  params.addRequiredParam<Real>("B", "'B' parameter for the psic func");
  params.addRequiredParam<Real>("C", "'C' parameter for the psic func");
  params.addRequiredParam<Real>("D", "'D' parameter");
  params.addRequiredParam<Real>("eta", "'eta' parameter");
  params.addRequiredParam<Real>("psic_orig", "'psic_oric'");
  params.addParam<MaterialPropertyName>(
      "fracture_toughness", "Gc", "The fracture toughness $\\Gc$");
  params.addParam<MaterialPropertyName>(
      "normalization_constant", "c0", "The normalization constant $c_0$");
  params.addParam<MaterialPropertyName>(
      "regularization_length", "l", "The phase-field regularization length");
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  return params;
}

PsicDeg::PsicDeg(const InputParameters & parameters)
  : PsicDegModel(parameters),
    _psie_active(getADMaterialProperty<Real>("psie_active")),
    _Gc(getADMaterialProperty<Real>(prependBaseName("fracture_toughness", true))),
    _c0(getADMaterialProperty<Real>(prependBaseName("normalization_constant", true))),
    _l(getADMaterialProperty<Real>(prependBaseName("regularization_length", true))),
    _eta(getParam<Real>("eta")),
    _A(getParam<Real>("A")),
    _B(getParam<Real>("B")),
    _C(getParam<Real>("C")),
    _D(getParam<Real>("D")),
    _d(coupledValue("phase_field")),
    _psic_orig(getParam<Real>("psic_orig"))
{
}
ADReal
PsicDeg::elasticEnergy(const ADReal & ep, const unsigned int derivative)
{

  ADReal result = 0;
  ADReal coef = _Gc[_qp] / _c0[_qp] / _l[_qp];
  // ADReal coef = _Gc;

  ADReal d = _d[_qp];

  ADReal dg1 = 1;
  ADReal dg2 = 2;
  ADReal dg = _psic_orig * ((_A * atan(_B * ep - _C) + _D));
  if (derivative == 1 || derivative == 2)
  {
    dg1 = _psic_orig * (_A * _B) / (std::pow(_C - _B * ep, 2) + 1);
  }
  if (derivative == 2)
  {
    dg2 = _psic_orig * (2 * _A * _B * _B * (_C - _B * ep)) /
          std::pow(std::pow(_C - _B * ep, 2) + 1, 2);
    // std::cout << "dg2 " << MetaPhysicL::raw_value(_C - _B * ep) << std::endl;
  }
  if (derivative == 1)
  {

    result = -(coef * (1 - d) * d * (_eta - 1) * dg1) / std::pow((coef * d - d * dg + dg), 2);
  }
  if (derivative == 2)
  {

    result = (coef * (d - 1) * d * (_eta - 1) *
              (coef * d * dg2 + (d - 1) * (2 * dg1 * dg1 - dg * dg2))) /
             std::pow(coef * d - (d - 1) * dg, 3);
  }
  return result * _psie_active[_qp];
}
