
#include "IrreversableNodalKernel.h"

registerMooseObject("MooseApp", IrreversableNodalKernel);

InputParameters
IrreversableNodalKernel::validParams()
{
  InputParameters params = NodalKernel::validParams();
  params.addClassDescription("Used to prevent a coupled variable from going below a lower bound");
  params.addRequiredCoupledVar(
      "v", "The coupled variable we require to be greater than the lower bound");
  params.addParam<std::vector<BoundaryName>>(
      "exclude_boundaries",
      {},
      "Boundaries on which not to execute the NodalKernel. This can be useful for avoiding "
      "singuarility in the matrix in case a constraint is active in the same place that a "
      "DirichletBC is set");
  return params;
}

IrreversableNodalKernel::IrreversableNodalKernel(const InputParameters & parameters)
  : NodalKernel(parameters),
    _v_var(coupled("v")),
    _v(coupledValue("v")),
    _v_old(coupledValueOld("v"))
{
  if (_var.number() == _v_var)
    mooseError("Coupled variable 'v' needs to be different from 'variable' with "
               "IrreversableNodalKernel");

  const auto & bnd_names = getParam<std::vector<BoundaryName>>("exclude_boundaries");
  for (const auto & bnd_name : bnd_names)
    _bnd_ids.insert(_mesh.getBoundaryID(bnd_name));
}

Real
IrreversableNodalKernel::computeQpResidual()
{
  for (auto bnd_id : _bnd_ids)
    if (_mesh.isBoundaryNode(_current_node->id(), bnd_id))
      return _u[_qp];

  return std::min(_u[_qp], _v[_qp] - _v_old[_qp]);
}

Real
IrreversableNodalKernel::computeQpJacobian()
{
  for (auto bnd_id : _bnd_ids)
    if (_mesh.isBoundaryNode(_current_node->id(), bnd_id))
      return 1;

  if (_u[_qp] <= _v[_qp] - _v_old[_qp])
    return 1;
  return 0;
}

Real
IrreversableNodalKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  for (auto bnd_id : _bnd_ids)
    if (_mesh.isBoundaryNode(_current_node->id(), bnd_id))
      return 0;

  if (jvar == _v_var)
    if (_v[_qp] - _v_old[_qp] < _u[_qp])
      return 1;

  return 0.0;
}
