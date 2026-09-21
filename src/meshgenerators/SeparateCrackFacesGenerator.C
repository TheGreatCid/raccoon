#include "SeparateCrackFacesGenerator.h"
#include "MooseMeshUtils.h"

#include "libmesh/boundary_info.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"

registerMooseObject("raccoonApp", SeparateCrackFacesGenerator);

InputParameters
SeparateCrackFacesGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundaries",
      "The crack face boundaries, e.g. the two sidesets BreakMeshByBlockGenerator creates. "
      "Boundaries that do not exist or are empty are skipped.");
  params.addRequiredRangeCheckedParam<Real>(
      "separation", "separation > 0", "Distance between the two crack faces after separation");
  params.addParam<bool>(
      "remove_boundaries",
      true,
      "Remove the crack face boundaries afterwards. A later cut then finds the crack from its "
      "geometry, and the faces of an earlier crack are not moved a second time.");
  // A plain name rather than a MeshGeneratorName: only the generator's meta data is read, the
  // generator must already be upstream of 'input'
  params.addParam<std::string>(
      "moved_nodes_from",
      "A mesh generator upstream of 'input' that moved nodes and recorded it (e.g. "
      "PlanarCrackGenerator); its offsets are combined with the ones made here");
  params.addClassDescription(
      "Moves the two faces of a broken crack apart by a small distance, leaving the crack front "
      "closed, so the crack is a thin slit rather than coincident faces.");
  return params;
}

SeparateCrackFacesGenerator::SeparateCrackFacesGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _boundaries(getParam<std::vector<BoundaryName>>("boundaries")),
    _separation(getParam<Real>("separation"))
{
  declareMeshProperty<std::vector<Point>>("moved_node_positions", {});
  declareMeshProperty<std::vector<Point>>("moved_node_offsets", {});
}

std::unique_ptr<MeshBase>
SeparateCrackFacesGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);
  if (!mesh->is_prepared())
    mesh->find_neighbors();
  auto & boundary_info = mesh->get_boundary_info();

  // Original positions recorded upstream, keyed by the nodes' current positions
  std::map<std::array<Real, 3>, Point> upstream;
  const auto key = [](const Point & p)
  {
    const auto round = [](const Real x) { return std::round(x * 1e10) / 1e10; };
    return std::array<Real, 3>{round(p(0)), round(p(1)), round(p(2))};
  };
  if (isParamValid("moved_nodes_from"))
  {
    const auto & generator = getParam<std::string>("moved_nodes_from");
    const auto & positions = getMeshProperty<std::vector<Point>>("moved_node_positions", generator);
    const auto & offsets = getMeshProperty<std::vector<Point>>("moved_node_offsets", generator);
    for (const auto i : index_range(positions))
      upstream[key(positions[i])] = positions[i] + offsets[i];
  }

  // Averaged outward normal of each node's crack faces, and the boundaries it lies on
  std::map<dof_id_type, RealVectorValue> normal_sum;
  std::map<dof_id_type, std::set<boundary_id_type>> node_boundaries;
  std::size_t n_faces = 0;
  std::set<boundary_id_type> crack_ids;
  for (const auto & name : _boundaries)
  {
    const auto id = MooseMeshUtils::getBoundaryID(name, *mesh);
    if (id != Moose::INVALID_BOUNDARY_ID)
      crack_ids.insert(id);
  }
  {
    for (const auto & [elem_id, side, id] : boundary_info.build_side_list())
    {
      if (!crack_ids.count(id))
        continue;
      const Elem * elem = mesh->elem_ptr(elem_id);
      const auto face = elem->build_side_ptr(side);
      if (face->n_vertices() < 3)
        continue;
      auto normal = (face->point(1) - face->point(0)).cross(face->point(2) - face->point(0));
      // Outward: pointing away from the element
      if (normal * (face->vertex_average() - elem->vertex_average()) < 0)
        normal = -normal;
      for (const auto n : elem->nodes_on_side(side))
      {
        normal_sum[elem->node_id(n)] += normal;
        node_boundaries[elem->node_id(n)].insert(id);
      }
      ++n_faces;
    }
  }

  std::map<dof_id_type, Point> original_positions;
  const auto original_of = [&](const Node & node)
  {
    const auto it = upstream.find(key(node));
    return it == upstream.end() ? Point(node) : it->second;
  };

  // Move the vertices of the crack faces; nodes on both faces (the crack front) stay
  std::set<dof_id_type> moved;
  std::size_t n_front = 0;
  for (const auto & [node_id, normal] : normal_sum)
  {
    if (node_boundaries[node_id].size() > 1)
    {
      ++n_front;
      continue;
    }
    if (normal.norm() == 0)
      continue;
    auto & node = mesh->node_ref(node_id);
    original_positions.emplace(node_id, original_of(node));
    node -= 0.5 * _separation * normal.unit();
    moved.insert(node_id);
  }

  // Straighten the edges of the elements around moved nodes, which also places the mid-edge nodes of
  // the crack faces, and check that none inverted
  std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elems;
  MeshTools::build_nodes_to_elem_map(*mesh, nodes_to_elems);
  std::set<dof_id_type> around;
  for (const auto node_id : moved)
    for (const auto * elem : nodes_to_elems[node_id])
      around.insert(elem->id());
  for (const auto elem_id : around)
  {
    Elem * elem = mesh->elem_ptr(elem_id);
    if (elem->default_order() == SECOND)
      for (const auto e : elem->edge_index_range())
      {
        const auto nodes = elem->nodes_on_edge(e);
        if (nodes.size() < 3)
          continue;
        auto & middle = elem->node_ref(nodes[2]);
        const Point midpoint = 0.5 * (elem->point(nodes[0]) + elem->point(nodes[1]));
        if ((midpoint - middle).norm() > 1e-12 * _separation)
        {
          original_positions.emplace(middle.id(), original_of(middle));
          middle = midpoint;
        }
      }
    // Jacobian of the straight-sided element at its vertices
    if (elem->dim() == 3 && elem->n_vertices() == 4)
    {
      const Point & a = elem->point(0);
      if ((elem->point(1) - a) * ((elem->point(2) - a).cross(elem->point(3) - a)) <= 0)
        mooseError("Separating the crack faces by ",
                   _separation,
                   " inverts element ",
                   elem->id(),
                   " at ",
                   elem->vertex_average(),
                   "; use a smaller separation.");
    }
  }

  auto & positions = setMeshProperty<std::vector<Point>>("moved_node_positions");
  auto & offsets = setMeshProperty<std::vector<Point>>("moved_node_offsets");
  positions.clear();
  offsets.clear();
  for (const auto & [node_id, original] : original_positions)
  {
    positions.push_back(mesh->point(node_id));
    offsets.push_back(original - mesh->point(node_id));
  }
  // Nodes moved upstream but not here keep their recorded offsets
  for (const auto & node : mesh->node_ptr_range())
    if (!original_positions.count(node->id()))
    {
      const auto it = upstream.find(key(*node));
      if (it != upstream.end())
      {
        positions.push_back(*node);
        offsets.push_back(it->second - *node);
      }
    }

  _console << name() << ": separated " << n_faces << " crack faces by " << _separation << ": moved "
           << moved.size() << " vertices, " << n_front << " crack front nodes left in place."
           << std::endl;
  if (n_faces == 0)
    mooseWarning("No faces on the crack boundaries; nothing was separated.");
  if (getParam<bool>("remove_boundaries"))
    for (const auto id : crack_ids)
      boundary_info.remove_id(id);

  mesh->unset_is_prepared();
  return mesh;
}
