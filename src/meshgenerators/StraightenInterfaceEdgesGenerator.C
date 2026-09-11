#include "StraightenInterfaceEdgesGenerator.h"
#include "MooseMeshUtils.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"
#include "libmesh/remote_elem.h"

registerMooseObject("raccoonApp", StraightenInterfaceEdgesGenerator);

InputParameters
StraightenInterfaceEdgesGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<std::vector<std::vector<SubdomainName>>>(
      "block_pairs",
      "The subdomain pairs whose shared faces are straightened, e.g. '10 11; 10 12'. Use the "
      "same pairs as the downstream BreakMeshByBlockGenerator.");
  params.addClassDescription(
      "Moves the mid-edge nodes of second-order faces shared by the listed block pairs to the "
      "straight-edge midpoint so BreakMeshByBlockGenerator can re-pair the broken sides of "
      "curved meshes. Nodes away from those faces are left untouched.");
  return params;
}

StraightenInterfaceEdgesGenerator::StraightenInterfaceEdgesGenerator(
    const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _block_pairs(getParam<std::vector<std::vector<SubdomainName>>>("block_pairs"))
{
  for (const auto & pair : _block_pairs)
    if (pair.size() != 2)
      paramError("block_pairs",
                 "Each entry must contain exactly two subdomains; got ",
                 pair.size(),
                 ".");
}

std::unique_ptr<MeshBase>
StraightenInterfaceEdgesGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  // Store each pair as (min, max) so the lookup is independent of which side we visit from.
  std::set<std::pair<subdomain_id_type, subdomain_id_type>> pairs;
  for (const auto & pair : _block_pairs)
  {
    const auto ids = MooseMeshUtils::getSubdomainIDs(*mesh, pair);
    for (const auto i : index_range(ids))
      if (!MooseMeshUtils::hasSubdomainID(*mesh, ids[i]))
        paramError("block_pairs", "The subdomain '", pair[i], "' was not found in the mesh.");
    pairs.emplace(std::min(ids[0], ids[1]), std::max(ids[0], ids[1]));
  }

  // The interface is found through neighbor links, which upstream generators may have dropped.
  if (!mesh->is_prepared())
    mesh->find_neighbors();

  for (const auto & elem : mesh->active_element_ptr_range())
    for (const auto s : elem->side_index_range())
    {
      const Elem * neighbor = elem->neighbor_ptr(s);
      if (!neighbor || neighbor == remote_elem)
        continue;

      const auto a = elem->subdomain_id();
      const auto b = neighbor->subdomain_id();
      if (!pairs.count({std::min(a, b), std::max(a, b)}))
        continue;

      // The side and its edges are proxies that share the mesh's Node objects, so moving a node
      // here moves it for every element around that edge and the mesh stays conforming.
      auto side = elem->build_side_ptr(s);
      if (side->n_nodes() > side->n_vertices() + side->n_edges())
        mooseError(name(),
                   ": interface face of element ",
                   elem->id(),
                   " has face-interior nodes (",
                   side->n_nodes(),
                   " nodes); only mid-edge nodes are straightened, so this element type "
                   "is not supported.");

      for (const auto e : side->edge_index_range())
      {
        auto edge = side->build_edge_ptr(e);
        // First-order edges have no mid-edge node to move.
        if (edge->n_nodes() < 3)
          continue;
        // EDGE3: nodes 0 and 1 are the vertices, node 2 is the mid-edge node.
        *edge->node_ptr(2) = 0.5 * (edge->point(0) + edge->point(1));
      }
    }

  mesh->unset_is_prepared();
  return mesh;
}
