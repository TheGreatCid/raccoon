//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPFFSourceRate.h"

registerMooseObject("raccoonApp", ADPFFSourceRate);

InputParameters
ADPFFSourceRate::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription("The source term in the phase-field evolution equation. The weak form "
                             "is $(w, \\dfrac{\\partial \\psi}{\\partial d})$.");
  params.addParam<MaterialPropertyName>("free_energy",
                                        "psi"
                                        "The Helmholtz free energy");
  params.addRequiredCoupledVar("phase_field_ddot",
                               "Name of the phase-field (damage) rate variable");

  return params;
}

ADPFFSourceRate::ADPFFSourceRate(const InputParameters & parameters)
  : ADKernelValue(parameters),
    DerivativeMaterialPropertyNameInterface(),
    _psi_name(getParam<MaterialPropertyName>("free_energy")),
    _ddot_name(getVar("phase_field_ddot", 0)->name()),
    _dpsi_dd_ddot(getADMaterialProperty<Real>(
        derivativePropertyNameSecond(_psi_name, _var.name(), {_ddot_name})))
{
}

ADReal
ADPFFSourceRate::precomputeQpResidual()
{
  return _dpsi_dd_ddot[_qp];
}
