#pragma once

#include "MeshGenerator.h"

/**
 * Removes the disjoint neighbor boundary pairs that BreakMeshByBlockGenerator registers, so the
 * two sides of a broken interface become ordinary, unconnected external boundaries.
 *
 * BreakMeshByBlockGenerator pairs each broken interface with its opposite side so elements keep
 * neighbor links across it (used by interface kernels such as cohesive zone models). libMesh
 * restores those links by locating the vertex average of each side on the opposite side, which
 * fails with "Periodic boundary neighbor not found" when the sides are curved, e.g. on a
 * second-order mesh moved into a deformed configuration. A traction-free crack needs no links
 * across it, so removing the pairs avoids the failure without altering the geometry.
 *
 * Do not use this when interface kernels or other objects need neighbors across the interface.
 */
class RemoveDisjointNeighborPairsGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  RemoveDisjointNeighborPairsGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// The mesh to modify
  std::unique_ptr<MeshBase> & _input;
};
