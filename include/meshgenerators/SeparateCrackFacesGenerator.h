#pragma once

#include "MeshGenerator.h"

/**
 * Pushes the two faces of a broken crack apart by a small distance, so tools that merge coincident
 * nodes (e.g. Cubit importing the mesh as geometry) see a thin slit instead of one surface.
 *
 * Every node on the given crack boundaries moves half the separation along the averaged outward
 * normal of its crack faces, reversed, i.e. away from the other face. Nodes on more than one of the
 * boundaries (the crack front, where the faces meet) stay put, so the slit closes at the front.
 * Mid-edge nodes of the elements around moved nodes are placed at their edge midpoints.
 *
 * The positions and offsets of moved nodes are stored as mesh meta data (moved_node_positions,
 * moved_node_offsets). With moved_nodes_from, the offsets recorded by an earlier generator (e.g.
 * PlanarCrackGenerator) are carried over and combined, so CopySolutionFields can look fields up at
 * the nodes' positions before either moved them.
 */
class SeparateCrackFacesGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  SeparateCrackFacesGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// The mesh to modify
  std::unique_ptr<MeshBase> & _input;
  /// Crack face boundaries (one per side)
  const std::vector<BoundaryName> & _boundaries;
  /// Distance between the two crack faces after separation
  const Real _separation;
};
