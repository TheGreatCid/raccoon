//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "NodeInterpolatedSolutionIC.h"

#include "libmesh/fe_interface.h"
#include "libmesh/fe_map.h"

registerMooseObject("raccoonApp", NodeInterpolatedSolutionIC);

InputParameters
NodeInterpolatedSolutionIC::validParams()
{
  InputParameters params = SolutionIC::validParams();
  params.addClassDescription(
      "Sets the initial condition from a field variable stored in an Exodus file, read only at the "
      "nodes of each element: values between the nodes are interpolated with the element's shape "
      "functions. Avoids looking up points on curved second-order elements, which can fail.");
  return params;
}

NodeInterpolatedSolutionIC::NodeInterpolatedSolutionIC(const InputParameters & parameters)
  : SolutionIC(parameters), _node_values_elem(nullptr)
{
}

bool
NodeInterpolatedSolutionIC::referencePoint(const Point & p, Point & reference) const
{
  // Three-dimensional elements only (see value())
  const unsigned int dim = 3;
  reference = _current_elem->reference_elem()->vertex_average();

  // Newton iteration on the element map. The map of a curved element is not injective beyond the
  // element, so starting at the centroid (rather than at a vertex) keeps the iteration on the
  // preimage inside the element.
  for (unsigned int it = 0; it < 20; ++it)
  {
    const Point residual = p - libMesh::FEMap::map(dim, _current_elem, reference);
    RealTensorValue J;
    for (const auto j : make_range(dim))
    {
      const Point dx = libMesh::FEMap::map_deriv(dim, _current_elem, j, reference);
      for (const auto i : make_range(dim))
        J(i, j) = dx(i);
    }
    Point step;
    J.solve(residual, step);
    reference += step;
    // Converged to well below the size of the reference element
    if (step.norm() < TOLERANCE * TOLERANCE)
      return libMesh::FEInterface::on_reference_element(reference, _current_elem->type(), TOLERANCE);
  }
  return false;
}

Real
NodeInterpolatedSolutionIC::value(const Point & p)
{
  // At a node, the lookup is the node itself
  if (_current_node || _current_elem->dim() != 3)
    return SolutionIC::value(p);

  Point reference;
  if (!referencePoint(p, reference))
    return SolutionIC::value(p);

  if (_node_values_elem != _current_elem)
  {
    _node_values.resize(_current_elem->n_nodes());
    for (const auto i : _current_elem->node_index_range())
      _node_values[i] = SolutionIC::value(_current_elem->point(i));
    _node_values_elem = _current_elem;
  }

  const libMesh::FEType fe_type(_current_elem->default_order(), libMesh::LAGRANGE);
  Real value = 0;
  for (const auto i : _current_elem->node_index_range())
    value += libMesh::FEInterface::shape(fe_type, _current_elem, i, reference) * _node_values[i];
  return value;
}
