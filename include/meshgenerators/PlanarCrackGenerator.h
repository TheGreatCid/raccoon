#pragma once

#include "MeshGenerator.h"

class PhaseFieldSolution;

/**
 * Turns a phase-field crack into a planar sharp crack: fits a plane to the ridge of the phase field,
 * moves nodes onto it so element faces tile the plane, and assigns the elements on the two sides of
 * the crack part of the plane to two subdomains for BreakMeshByBlockGenerator to split.
 *
 * 1. Band: elements whose centroid value reaches the threshold.
 * 2. Plane: the normal is given or taken from the band's structure tensor sum(V grad(d) grad(d)^T);
 *    the ridge (largest d along the normal) is found through every band centroid, and the plane
 *    passes through the mean of those ridge points weighted by their peak value.
 * 3. Crack extent: a point on the plane is inside the crack when the ridge there reaches the
 *    threshold.
 * 4. Snapping: every edge of an element inside the extent that crosses the plane gets one end moved
 *    onto the plane, the nearer one unless that would degrade an element below the quality limit.
 *    If neither end can move directly, free nodes around it are smoothed to make room. Nodes on
 *    external boundaries slide within them. Mid-edge nodes of the affected elements go to
 *    their edge midpoints, so the crack faces are exactly planar.
 * 5. Sides: the elements of faces lying on the plane inside the extent, and elements there that
 *    still cross the plane because their edges could not be snapped, go to the subdomain of the side
 *    their centroid is on; around unsnapped elements the crack follows element faces (sawtooth)
 *    instead of leaving a hole. Every node between the two sides where the ridge reaches the
 *    threshold then takes all elements around it; nodes where it does not stay shared, forming the
 *    crack front.
 *
 * The positions and offsets of moved nodes are stored as mesh meta data (moved_node_positions,
 * moved_node_offsets) so CopySolutionFields can look fields up at the original positions.
 * Tetrahedral meshes only; crack branching and curved cracks are not represented.
 */
class PlanarCrackGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  PlanarCrackGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  /// Largest ridge value over a point's projection onto the plane and the given extra points'
  Real crackRidgeValue(PhaseFieldSolution & phase_field,
                       const std::vector<Point> & points,
                       Real reach) const;

  /// The mesh to modify
  std::unique_ptr<MeshBase> & _input;
  /// Elements whose centroid value reaches this form the band; the crack covers the part of the
  /// plane where the ridge reaches it
  const Real _threshold;
  /// Search distance for the ridge on each side of a point (0: twice the largest band element)
  const Real _ridge_search_distance;
  /// Moves may not bring the worst surrounding scaled Jacobian below this (or the prior worst)
  const Real _min_scaled_jacobian;
  /// Smallest cosine between a boundary node's slide direction and the crack normal
  const Real _min_boundary_alignment;
  /// Element layers of free nodes smoothed when a node cannot be snapped directly (0: none)
  const unsigned int _relaxation_layers;
  /// Smoothing sweeps per relaxed move
  const unsigned int _relaxation_iterations;
  /// Subdomain for elements on the side opposite to the normal
  const subdomain_id_type _lower_block_id;
  /// Subdomain for elements on the side the normal points to
  const subdomain_id_type _upper_block_id;

  /// Crack plane point and unit normal, set by generate()
  Point _plane_point;
  RealVectorValue _plane_normal;
};
