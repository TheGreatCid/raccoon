//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPFFDiffusion.h"
#include "ADRankTwoTensorForward.h"
#include "Adaptivity.h"
#include "Assembly.h"
#include "EigenADReal.h"
#include "MooseTypes.h"
#include "MooseVariableField.h"
#include "SolutionUserObject.h"
#include "AuxiliarySystem.h"
#include "libmesh/id_types.h"

registerMooseObject("raccoonApp", ADPFFDiffusion);

InputParameters
ADPFFDiffusion::validParams()
{
  InputParameters params = ADKernel::validParams();
  params += BaseNameInterface::validParams();
  params.addClassDescription("The diffusion term in the phase-field evolution equation. The weak "
                             "form is $(\\grad w, \\dfrac{2\\Gc l}{c_0} \\grad d)$.");

  params.addParam<MaterialPropertyName>(
      "fracture_toughness", "Gc", "The fracture toughness $\\Gc$");
  params.addParam<MaterialPropertyName>(
      "normalization_constant", "c0", "The normalization constant $c_0$");
  params.addParam<MaterialPropertyName>(
      "regularization_length", "l", "The phase-field regularization length");
  params.addParam<bool>("recover", false, "Are you trying to recover");
  params.addCoupledVar("d_old", "d from recover");
  params.addCoupledVar("d_old_grad_ref", "d from recover");
  params.addCoupledVar("grad_dx", "d from recover");
  params.addCoupledVar("grad_dy", "d from recover");
  params.addCoupledVar("grad_dz", "d from recover");
  params.addCoupledVar("Fnobar_xx_1", "F_init from recover");
  params.addCoupledVar("Fnobar_xx_2", "F_init from recover");

  params.addCoupledVar("F_current", "Current F from recover");

  return params;
}

ADPFFDiffusion::ADPFFDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    BaseNameInterface(parameters),
    _Gc(getADMaterialProperty<Real>(prependBaseName("fracture_toughness", true))),
    _c0(getADMaterialProperty<Real>(prependBaseName("normalization_constant", true))),
    _l(getADMaterialProperty<Real>(prependBaseName("regularization_length", true))),
    _recover(getParam<bool>("recover")),
    _d_old_grad(coupledGradient("d_old")),
    _d_old_grad_ref(coupledValue("d_old_grad_ref")),
    _grad_dx(coupledValue("grad_dx")),
    _grad_dy(coupledValue("grad_dy")),
    _grad_dz(coupledValue("grad_dz")),
    _Fnobar_xx_1(coupledValue("Fnobar_xx_1")),
    _Fnobar_xx_2(coupledValue("Fnobar_xx_2")),
    _F_current(coupledValue("F_current"))

// _d_diff(dynamic_cast<MooseVariable &>(_subproblem.getVariable(0, "diff")))
// _d_diff(_sys.addVector("d_diff", true, GHOSTED))
{
}

ADReal
ADPFFDiffusion::computeQpResidual()
{
  // if (_sys.subproblem().nLinearIterations(_sys.number()) == 1)
  //   std::cout << _sys.subproblem().nLinearIterations(_sys.number()) << std::endl;
  if (!_recover)
  {
    ADReal value = _grad_test[_i][_qp] * _grad_u[_qp];
    // auto dofs = _var.dofIndices();

    // auto currentdof = dofs[_i];
    // _d_diff.set(currentdof, MetaPhysicL::raw_value(value));
    // auto currnode = _current_elem->node_ptr(_qp);

    // auto dofnum = currnode->dof_number(_sys.number(), 0, 0);
    // _d_diff.setNodalValue(MetaPhysicL::raw_value(value));
    // vec->set(dofnum, MetaPhysicL::raw_value(value));
    // std::cout << " =====" << std::endl;

    return 2 * _Gc[_qp] * _l[_qp] / _c0[_qp] * value;
  }
  else
  {
    ADVariableGradient gra(3);

    gra.setAllValues({_grad_dx[_qp], 0, 0});

    auto dot = _grad_test[_i][_qp](0) * _grad_dx[_qp]; // + _grad_test[_i][_qp](1) * _grad_dy[_qp] +
                                                       //  _grad_test[_i][_qp](2) * _grad_dz[_qp];

    ADReal value = 0;

    if (_qp == 0)
      value = _F_current[_qp] * _Fnobar_xx_1[_qp];
    if (_qp == 1)
      value = _F_current[_qp] * _Fnobar_xx_2[_qp];

    value *= _grad_test[_i][_qp] * _grad_u[_qp];
    // ADReal value = dot + _grad_test[_i][_qp] * (_grad_u[_qp] - _d_old_grad[_qp]);
    return 2 * _Gc[_qp] * _l[_qp] / _c0[_qp] * value;

    // MooseVariable & min = (_u[_qp] - _d_old[_qp]);
    // min.
  }
}
