#pragma once

#include "MeshGenerator.h"

/**
 * Renames node sets independently of their side sets.
 *
 * libMesh stores side-set and node-set names in separate maps, even when a boundary's
 * side set and node set share the same boundary id. MOOSE normally derives node sets from
 * side sets (construct_node_list_from_side_list) and gives them the same name, which some
 * downstream tools cannot disambiguate. This generator changes only the node-set name map,
 * leaving the side-set names untouched, so a boundary's node set and side set end up with
 * different names.
 *
 * Set the [Mesh] block's `construct_node_list_from_side_list = false` so the rename applied
 * here is not overwritten when the mesh is finalized.
 */
class RenameNodesetGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  RenameNodesetGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// The mesh to modify
  std::unique_ptr<MeshBase> & _input;
  /// Boundaries whose node-set name will be changed (side-set names are kept)
  const std::vector<BoundaryName> _old_nodesets;
  /// New node-set names, one per entry in _old_nodesets
  const std::vector<BoundaryName> _new_nodesets;
  /// Build node lists from the side lists before renaming
  const bool _build_node_list;
};
