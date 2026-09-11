#pragma once

#include "MeshGenerator.h"

/**
 * Moves the mid-edge nodes of second-order faces shared by the listed block pairs to the
 * straight-edge midpoint, leaving every other node untouched.
 *
 * BreakMeshByBlockGenerator registers the two sides of each broken interface as a disjoint
 * (zero-translation periodic) boundary pair. libMesh re-pairs those sides by locating the
 * vertex average of each side on the opposite side. On a curved second-order face that point
 * lies slightly off the face, the opposite element is not found, and mesh preparation fails
 * with "Periodic boundary neighbor not found". Straightening only the interface faces before
 * breaking avoids that while keeping curved geometry (e.g. circular boundaries) elsewhere.
 *
 * Place this generator between the subdomain assignment and BreakMeshByBlockGenerator, with
 * the same block_pairs.
 */
class StraightenInterfaceEdgesGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  StraightenInterfaceEdgesGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// The mesh to modify
  std::unique_ptr<MeshBase> & _input;
  /// Block pairs whose shared faces are straightened
  const std::vector<std::vector<SubdomainName>> _block_pairs;
};
