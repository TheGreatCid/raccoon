#include "RenameNodesetGenerator.h"
#include "MooseMeshUtils.h"

#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"

registerMooseObject("raccoonApp", RenameNodesetGenerator);

InputParameters
RenameNodesetGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "old_nodeset",
      "The boundaries whose node-set name will be changed. Their side-set names are left "
      "untouched, so the node set and side set end up with different names.");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "new_nodeset", "The new node-set names, one per entry in 'old_nodeset'.");
  params.addParam<bool>(
      "build_node_list_from_side_list",
      true,
      "Build node lists from the side lists first so node sets exist for side-set-only "
      "boundaries. Set the [Mesh] block's 'construct_node_list_from_side_list = false' so the "
      "rename performed here is not overwritten when the mesh is finalized.");
  params.addClassDescription(
      "Renames node sets independently of their side sets so a boundary's node-set name can "
      "differ from its side-set name (libMesh stores them in separate name maps).");
  return params;
}

RenameNodesetGenerator::RenameNodesetGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _old_nodesets(getParam<std::vector<BoundaryName>>("old_nodeset")),
    _new_nodesets(getParam<std::vector<BoundaryName>>("new_nodeset")),
    _build_node_list(getParam<bool>("build_node_list_from_side_list"))
{
  if (_old_nodesets.size() != _new_nodesets.size())
    paramError("new_nodeset",
               "The number of 'new_nodeset' names (",
               _new_nodesets.size(),
               ") must match the number of 'old_nodeset' names (",
               _old_nodesets.size(),
               ").");
}

std::unique_ptr<MeshBase>
RenameNodesetGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);
  auto & boundary_info = mesh->get_boundary_info();

  // Make sure node sets exist for side-set-only boundaries before renaming.
  if (_build_node_list)
    boundary_info.build_node_list_from_side_list();

  for (std::size_t i = 0; i < _old_nodesets.size(); ++i)
  {
    const auto id = MooseMeshUtils::getBoundaryID(_old_nodesets[i], *mesh);
    if (id == libMesh::BoundaryInfo::invalid_id)
      paramError("old_nodeset", "The boundary '", _old_nodesets[i], "' was not found in the mesh.");

    // Only the node-set name map is touched; the side-set name for this id is left as-is.
    boundary_info.nodeset_name(id) = _new_nodesets[i];
  }

  mesh->set_isnt_prepared();
  return mesh;
}
