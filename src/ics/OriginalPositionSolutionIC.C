#include "OriginalPositionSolutionIC.h"

#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_map.h"

registerMooseObject("raccoonApp", OriginalPositionSolutionIC);

InputParameters
OriginalPositionSolutionIC::validParams()
{
  InputParameters params = SolutionIC::validParams();
  params.addClassDescription("Sets the initial condition from a field variable stored in an Exodus "
                             "file, looked up at the positions nodes had before they were moved.");
  params.addRequiredParam<std::vector<Point>>("moved_node_positions",
                                              "Current positions of the moved nodes");
  params.addRequiredParam<std::vector<Point>>(
      "moved_node_offsets", "Original minus current position of each moved node");
  return params;
}

OriginalPositionSolutionIC::OriginalPositionSolutionIC(const InputParameters & parameters)
  : SolutionIC(parameters),
    _moved_positions(getParam<std::vector<Point>>("moved_node_positions")),
    _moved_offsets(getParam<std::vector<Point>>("moved_node_offsets"))
{
  if (_moved_positions.size() != _moved_offsets.size())
    paramError("moved_node_offsets",
               "There must be one offset per position (",
               _moved_positions.size(),
               ").");
  if (!_moved_positions.empty())
    _tree = std::make_unique<KDTree>(_moved_positions, 10);
}

Point
OriginalPositionSolutionIC::offset(const Point & p)
{
  if (!_tree)
    return Point();
  std::vector<std::size_t> index;
  std::vector<Real> distance_sqr(1);
  _tree->neighborSearch(p, 1, index, distance_sqr);
  // Moved positions are copied exactly into the mesh, including onto nodes duplicated by a break
  return distance_sqr[0] <= Utility::pow<2>(TOLERANCE * TOLERANCE * std::max(1., p.norm()))
             ? _moved_offsets[index[0]]
             : Point();
}

Real
OriginalPositionSolutionIC::value(const Point & p)
{
  if (_current_node)
    return SolutionIC::value(p + offset(p));

  // A quadrature point: interpolate the element's nodal offsets at its reference coordinates
  const Point reference = libMesh::FEMap::inverse_map(_current_elem->dim(), _current_elem, p);
  const libMesh::FEType fe_type(_current_elem->default_order(), LAGRANGE);
  Point shift;
  for (const auto i : _current_elem->node_index_range())
    shift += libMesh::FEInterface::shape(fe_type, _current_elem, i, reference) *
             offset(_current_elem->point(i));
  return SolutionIC::value(p + shift);
}
