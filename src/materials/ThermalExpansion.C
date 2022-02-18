//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ThermalExpansion.h"

registerADMooseObject("raccoonApp", ThermalExpansion);

InputParameters
ThermalExpansion::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription(
      "This class computes the thermal deformation gradient, $\\bfF = (1+\\alpha(T-T_0))\bfI$, "
      "where \\alpha is the thermal expansion coefficient, $T$ is the temperature, and $T_0$ is "
      "the reference temperature corresponding to zero thermal expansion.");

  params.addRequiredParam<MaterialPropertyName>("eigen_deformation_gradient_name",
                                                "The name of the thermal deformation gradient");
  params.addRequiredParam<Real>("thermal_expansion_function",
                                "Function describing the thermal expansion coefficient $\alpha$");
  params.addRequiredCoupledVar("temperature", "The current temperature");
  params.addRequiredCoupledVar("reference_temperature",
                               "The reference temperature corresponding to zero thermal expansion");

  params.suppressParameter<bool>("use_displaced_mesh");
  return params;
}

ThermalExpansion::ThermalExpansion(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _Fg(declareADProperty<RankTwoTensor>(prependBaseName("eigen_deformation_gradient_name", true))),
    _Fg_old(getMaterialPropertyOldByName<RankTwoTensor>(
        prependBaseName("eigen_deformation_gradient_name", true))),
    _alpha(getParam<Real>("thermal_expansion_function")),
    _T(adCoupledValue("temperature")),
    _T_old(coupledValueOld("temperature")),
    _T0(coupledValue("reference_temperature"))
{
}

void
ThermalExpansion::initQpStatefulProperties()
{
  _Fg[_qp].setToIdentity();
}

void
ThermalExpansion::computeQpProperties()
{
  const Real old_temp = _t_step <= 1 ? _T0[_qp] : _T_old[_qp];

  // const Point p;
  // const Real alpha = _alpha.value(MetaPhysicL::raw_value(_T[_qp]), p);
  // const Real alpha_old = _alpha.value(old_temp, p);
  // const Real alpha_avg = (alpha + alpha_old) / 2;
  if (_dt != 0)
    _Fg[_qp].addIa(_alpha * (_T[_qp] - old_temp) / _dt);
}
