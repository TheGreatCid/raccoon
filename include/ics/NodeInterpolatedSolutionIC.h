//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "SolutionIC.h"

/**
 * A SolutionIC that reads the solution only at the nodes of each element.
 *
 * The initial condition of a second-order variable is sampled at the nodes and at quadrature
 * points along each element edge, from which the mid-edge degrees of freedom are projected. Looking
 * up an edge point in the solution requires inverting the curved element's map, and libMesh starts
 * that Newton iteration at a reference vertex: on some curved TET10 elements it converges to a
 * second preimage outside the element, and the point is reported as outside the mesh. Here an edge
 * point is instead evaluated by interpolating the solution at the current element's nodes, with the
 * point's reference coordinates found by a Newton iteration started at the element centroid. When
 * the mesh is the solution file's mesh (a restart from a remapped file), the solution on each
 * element is the same polynomial, so the interpolation reproduces it exactly.
 */
class NodeInterpolatedSolutionIC : public SolutionIC
{
public:
  static InputParameters validParams();

  NodeInterpolatedSolutionIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  /// Reference coordinates of a point in the current element, from a Newton iteration started at the
  /// reference centroid; false if it does not converge inside the element
  bool referencePoint(const Point & p, Point & reference) const;

  /// The solution at the nodes of the element they were last read for
  std::vector<Real> _node_values;
  /// The element _node_values belong to
  const Elem * _node_values_elem;
};
