#pragma once

#include "MeshGenerator.h"

/**
 * Splits the damaged band of a phase-field crack into two subdomains, one on each side of the
 * crack mid-surface, so that BreakMeshByBlockGenerator can open a sharp crack along it.
 *
 * The phase-field variable is read from an Exodus file and evaluated at every element centroid.
 * Elements with d >= threshold form the band. Within the band, grad(d) points toward the ridge
 * (the crack mid-surface), so the sign of grad(d) . n, with n the crack normal, tells which side
 * of the ridge an element lies on. The normal is either given, or estimated per element from the
 * structure tensor sum(grad(d) grad(d)^T) of nearby band elements and then oriented consistently
 * by walking across the band. Elements where grad(d) . n is too small to trust (e.g. a d = 1
 * plateau) are filled in from their neighbors, advancing from both flanks so the two sides meet
 * in the middle, and a few majority-vote passes remove isolated misassignments.
 *
 * Elements below the threshold keep their subdomain, so the break between the two new blocks
 * ends where the band ends and that end becomes the crack front. Crack branching is not
 * supported: a two-way split cannot represent a junction.
 *
 * Typical use: this generator, then StraightenInterfaceEdgesGenerator (for curved second-order
 * meshes), then BreakMeshByBlockGenerator, the last two with block_pairs = 'lower upper'.
 */
class PhaseFieldCrackSubdomainGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  PhaseFieldCrackSubdomainGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// The mesh to modify
  std::unique_ptr<MeshBase> & _input;
  /// Exodus file holding the phase-field variable
  const FileName & _solution_file;
  /// Name of the nodal phase-field variable in the file
  const std::string & _variable;
  /// Elements whose centroid value is at least this are part of the crack band
  const Real _threshold;
  /// Whether a fixed crack normal was given
  const bool _has_normal;
  /// Radius of the neighborhood used to estimate the crack normal
  const Real _normal_radius;
  /// |grad(d) . n| below this fraction of the band's largest |grad(d)| is treated as undecided
  const Real _gradient_tolerance;
  /// Number of majority-vote smoothing passes
  const unsigned int _smoothing_passes;
  /// Subdomain for band elements on the side of the ridge opposite to the normal
  const subdomain_id_type _lower_block_id;
  /// Subdomain for band elements on the side of the ridge the normal points to
  const subdomain_id_type _upper_block_id;
};
