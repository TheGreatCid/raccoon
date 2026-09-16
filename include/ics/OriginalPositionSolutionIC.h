#pragma once

#include "SolutionIC.h"
#include "KDTree.h"

/**
 * A SolutionIC for a mesh whose nodes were moved after the solution was written (e.g. by
 * PlanarCrackGenerator): values are looked up at each point's position before the move. Nodes are
 * shifted by their recorded offset; quadrature points by the element's nodal offsets interpolated
 * with its shape functions. Nodes without a recorded offset are taken as unmoved.
 */
class OriginalPositionSolutionIC : public SolutionIC
{
public:
  static InputParameters validParams();

  OriginalPositionSolutionIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  /// The recorded offset (original - current) of a node at position p, zero if none
  Point offset(const Point & p);

  /// Current positions of the moved nodes
  const std::vector<Point> & _moved_positions;
  /// Original minus current position of each moved node
  const std::vector<Point> & _moved_offsets;
  /// Search tree over the moved positions
  std::unique_ptr<KDTree> _tree;
};
