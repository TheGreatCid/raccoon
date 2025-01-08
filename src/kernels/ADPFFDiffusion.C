//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADPFFDiffusion.h"
#include "ADRankTwoTensorForward.h"
#include "Adaptivity.h"
#include "Assembly.h"
#include "EigenADReal.h"
#include "MooseTypes.h"
#include "MooseVariableFE.h"
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
  params.addParam<Real>("qp", 8, "Number of qps");
  params.addParam<UserObjectName>("solution", "The SolutionUserObject to extract data from.");
  // params.addCoupledVar("grad_xx", "xx gradient");
  // params.addCoupledVar("grad_xy", "xy gradient");
  // params.addCoupledVar("grad_xz", "xz gradient");
  // params.addCoupledVar("grad_yx", "xx gradient");
  // params.addCoupledVar("grad_yy", "xy gradient");
  // params.addCoupledVar("grad_yz", "xz gradient");
  // params.addCoupledVar("grad_zx", "xx gradient");
  // params.addCoupledVar("grad_zy", "xy gradient");
  // params.addCoupledVar("grad_zz", "xz gradient");
  return params;
}

ADPFFDiffusion::ADPFFDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    BaseNameInterface(parameters),
    _Gc(getADMaterialProperty<Real>(prependBaseName("fracture_toughness", true))),
    _c0(getADMaterialProperty<Real>(prependBaseName("normalization_constant", true))),
    _l(getADMaterialProperty<Real>(prependBaseName("regularization_length", true))),
    _recover(getParam<bool>("recover")),
    // _grad_disp(adCoupledGradients("displacements")),
    // _grad_xx(coupledValue("grad_xx")),
    // _grad_xy(coupledValue("grad_xy")),
    // _grad_xz(coupledValue("grad_xz")),
    // _grad_yx(coupledValue("grad_yx")),
    // _grad_yy(coupledValue("grad_yy")),
    // _grad_yz(coupledValue("grad_yz")),
    // _grad_zx(coupledValue("grad_zx")),
    // _grad_zy(coupledValue("grad_zy")),
    // _grad_zz(coupledValue("grad_zz")),
    _solution_object_ptr(NULL)

// _d_diff(dynamic_cast<MooseVariable &>(_subproblem.getVariable(0, "diff")))
// _d_diff(_sys.addVector("d_diff", true, GHOSTED))
{
}

void
ADPFFDiffusion::initialSetup()
{
  if (_recover == true)
    _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");

  // auto dim = _mesh.dimension();
  // for (unsigned i = dim; i < 3; ++i)
  // {
  //   _grad_disp.push_back(&_ad_grad_zero);
  // }
}

ADReal
ADPFFDiffusion::computeQpResidual()
{
  // If we are not recovering
  if (!_recover)
  {
    ADReal value = _grad_test[_i][_qp] * _grad_u[_qp];

    return 2 * _Gc[_qp] * _l[_qp] / _c0[_qp] * value;
  }
  // If we are recovering
  else
  {
    //======================================================
    // Calculating the current def grad
    // ADRankTwoTensor A;
    // A(0, 0) = _grad_xx[_qp];
    // A(0, 1) = _grad_xy[_qp];
    // A(0, 2) = _grad_xz[_qp];

    // A(1, 0) = _grad_yx[_qp];
    // A(1, 1) = _grad_yy[_qp];
    // A(1, 2) = _grad_yz[_qp];

    // A(2, 0) = _grad_zx[_qp];
    // A(2, 1) = _grad_zy[_qp];
    // A(2, 2) = _grad_zz[_qp];
    // auto F = A;
    // F.addIa(1.0);

    ADRankTwoTensor F_prev;
    F_prev.setToIdentity();
    // Getting past def grad
    std::vector<std::string> indices = {"x", "y", "z"};
    auto dim = _mesh.dimension();
    // Populate tensor from solution object
    for (unsigned int i_ind = 0; i_ind < dim; i_ind++)
      for (unsigned int j_ind = 0; j_ind < dim; j_ind++)
      {
        F_prev(i_ind, j_ind) = _solution_object_ptr->directValue(
            _current_elem,
            "Fnobar_" + indices[i_ind] + indices[j_ind] + "_" + std::to_string(_qp + 1));
      }
    // multipying in old def grad
    // F = F * F_prev;

    //======================================================

    // THE TEST GRADIENT AND _grad_disp NEED TO BE ON DIFFERENT CONFIGURATIONS

    ADReal value = 0;

    auto mult = (_grad_u[_qp] * F_prev);
    value = _grad_test[_i][_qp] * F_prev * mult;
    // ADReal value = dot + _grad_test[_i][_qp] * (_grad_u[_qp] - _d_old_grad[_qp]);
    return 2 * _Gc[_qp] * _l[_qp] / _c0[_qp] * value;

    // MooseVariable & min = (_u[_qp] - _d_old[_qp]);
    // min.
  }
}
