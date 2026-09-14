#include "RemoveDisjointNeighborPairsGenerator.h"
#include "MooseMeshUtils.h"

#include "libmesh/mesh_base.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary_base.h"

registerMooseObject("raccoonApp", RemoveDisjointNeighborPairsGenerator);

InputParameters
RemoveDisjointNeighborPairsGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addParam<std::vector<BoundaryName>>(
      "boundaries",
      "Only remove pairs involving these boundaries. By default every pair is removed.");
  params.addClassDescription(
      "Removes the disjoint neighbor boundary pairs registered by BreakMeshByBlockGenerator so the "
      "sides of a broken interface are unconnected external boundaries (e.g. a traction-free "
      "crack). Needed when those sides are curved, where libMesh cannot re-pair them.");
  return params;
}

RemoveDisjointNeighborPairsGenerator::RemoveDisjointNeighborPairsGenerator(
    const InputParameters & parameters)
  : MeshGenerator(parameters), _input(getMesh("input"))
{
}

std::unique_ptr<MeshBase>
RemoveDisjointNeighborPairsGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  std::set<boundary_id_type> restricted;
  if (isParamValid("boundaries"))
    for (const auto & boundary : getParam<std::vector<BoundaryName>>("boundaries"))
    {
      const auto id = MooseMeshUtils::getBoundaryID(boundary, *mesh);
      if (id == Moose::INVALID_BOUNDARY_ID)
        paramError("boundaries", "The boundary '", boundary, "' was not found in the mesh.");
      restricted.insert(id);
    }

  const auto * pairs = mesh->get_disjoint_neighbor_boundary_pairs();
  if (!pairs || pairs->empty())
  {
    mooseWarning(name(), ": the mesh has no disjoint neighbor boundary pairs to remove.");
    return mesh;
  }

  // Collect first: removing invalidates the map iterators.
  std::vector<std::pair<boundary_id_type, boundary_id_type>> to_remove;
  for (const auto & [id, boundary] : *pairs)
    if (restricted.empty() || restricted.count(boundary->myboundary) ||
        restricted.count(boundary->pairedboundary))
      to_remove.emplace_back(boundary->myboundary, boundary->pairedboundary);

  for (const auto & [b1, b2] : to_remove)
    mesh->remove_disjoint_boundary_pair(b1, b2);

  // Neighbor links across the removed pairs must be rebuilt without them.
  mesh->unset_has_neighbor_ptrs();
  return mesh;
}
