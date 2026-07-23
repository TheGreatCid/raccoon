//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NodalConstantBounds.h"
#include "SystemBase.h"
#include "PetscSupport.h"
#include "FEProblemBase.h"

registerMooseObject("raccoonApp", NodalConstantBounds);

InputParameters
NodalConstantBounds::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Provides a constant bound of a nodal variable for PETSc's variational inequalities "
      "solver. The dummy aux variable need not match the family of the bounded variable, "
      "allowing families like BERNSTEIN that the AuxKernel system rejects at second order.");
  MooseEnum type_options("upper=0 lower=1", "upper");
  params.addParam<MooseEnum>(
      "bound_type",
      type_options,
      "Type of bound. 'upper' refers to the upper bound. 'lower' refers to the lower value.");
  params.addRequiredParam<NonlinearVariableName>("bounded_variable", "The variable to be bounded");
  params.addRequiredParam<Real>("bound_value", "The value of bound for the variable");
  params.registerBase("Bounds");
  return params;
}

NodalConstantBounds::NodalConstantBounds(const InputParameters & parameters)
  : AuxKernel(parameters),
    _type((BoundType)(int)parameters.get<MooseEnum>("bound_type")),
    _bounded_vector(_type == UPPER ? _nl_sys.getVector("upper_bound")
                                   : _nl_sys.getVector("lower_bound")),
    _bounded_var(_nl_sys.getVariable(_tid, getParam<NonlinearVariableName>("bounded_variable"))),
    _bound_value(getParam<Real>("bound_value"))
{
  // The dummy must be nodal so we iterate over nodes, and of the same order as the
  // bounded variable so every node carrying a bounded-variable DoF is visited.
  if (!isNodal())
    paramError("variable", "The dummy bounds aux variable must be a nodal (Lagrange) variable.");
  if (!_bounded_var.isNodal())
    paramError("bounded_variable",
               "The bounded variable must have node-associated degrees of freedom.");
  if (_bounded_var.feType().order != _var.feType().order)
    paramError("variable",
               "The dummy bounds aux variable must have the same order as the bounded variable.");
}

Real
NodalConstantBounds::computeValue()
{
  const auto sys_num = _nl_sys.number();
  const auto var_num = _bounded_var.number();

  if (_current_node->n_dofs(sys_num, var_num) > 0)
  {
    mooseAssert(_current_node->n_dofs(sys_num, var_num) == 1,
                "NodalConstantBounds requires exactly one DoF per node for the bounded variable");
    const dof_id_type dof = _current_node->dof_number(sys_num, var_num, 0);
    _bounded_vector.set(dof, _bound_value);
  }

  return 0.0;
}
