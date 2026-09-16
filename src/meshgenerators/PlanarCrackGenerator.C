#include "PlanarCrackGenerator.h"
#include "MooseMeshUtils.h"
#include "PhaseFieldSolution.h"
#include "RankTwoTensor.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/remote_elem.h"
#include "libmesh/string_to_enum.h"

#include <optional>

registerMooseObject("raccoonApp", PlanarCrackGenerator);

namespace
{
/// Tetrahedron vertex pairs forming its edges
const std::array<std::array<unsigned int, 2>, 6> tet_edges = {
    {{{0, 1}}, {{1, 2}}, {{0, 2}}, {{0, 3}}, {{1, 3}}, {{2, 3}}}};

/**
 * Smallest corner scaled Jacobian of the straight-sided tetrahedron on the element's vertices,
 * normalized so a regular tetrahedron gives 1 (negative when inverted)
 */
Real
scaledJacobian(const Elem & elem)
{
  const Point & a = elem.point(0);
  const Point & b = elem.point(1);
  const Point & c = elem.point(2);
  const Point & d = elem.point(3);
  const Real six_volume = (b - a) * ((c - a).cross(d - a));
  Real largest_product = 0;
  for (unsigned int corner = 0; corner < 4; ++corner)
  {
    Real product = 1;
    for (const auto & edge : tet_edges)
      if (edge[0] == corner || edge[1] == corner)
        product *= (elem.point(edge[0]) - elem.point(edge[1])).norm();
    largest_product = std::max(largest_product, product);
  }
  return largest_product > 0 ? std::sqrt(2.) * six_volume / largest_product : 0;
}
}

InputParameters
PlanarCrackGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params += PhaseFieldSolution::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredRangeCheckedParam<Real>(
      "threshold",
      "threshold > 0 & threshold <= 1",
      "Elements whose phase-field value at the centroid reaches this form the band used to fit the "
      "plane; the crack covers the part of the plane where the ridge reaches it");
  params.addParam<RealVectorValue>(
      "normal",
      "A fixed crack plane normal; only the plane's position is fitted. The side it points to "
      "becomes 'upper_block_id'. By default the normal is fitted as well.");
  params.addRangeCheckedParam<Real>(
      "ridge_search_distance",
      0,
      "ridge_search_distance >= 0",
      "Distance on each side of a point searched along the normal for the ridge (largest d). 0 uses "
      "twice the largest band element size.");
  params.addRangeCheckedParam<Real>(
      "min_scaled_jacobian",
      0.2,
      "min_scaled_jacobian > 0 & min_scaled_jacobian < 1",
      "A node is not moved onto the plane if that brings the worst scaled Jacobian (1 for a regular "
      "tetrahedron) of the surrounding elements below this, or below the worst they already had if "
      "that was lower.");
  params.addRangeCheckedParam<Real>(
      "min_boundary_alignment",
      0.5,
      "min_boundary_alignment > 0 & min_boundary_alignment <= 1",
      "Nodes on external boundaries slide within the boundary onto the plane (along the face, or "
      "along the line where two faces meet). The slide is refused when the cosine between its "
      "direction and the normal is below this, which limits it to 1/alignment times the distance "
      "to the plane.");
  params.addParam<unsigned int>(
      "relaxation_layers",
      2,
      "When a node cannot be moved onto the plane without degrading an element, the nodes not on the "
      "plane within this many element layers are smoothed to make room (nodes on external "
      "boundaries only within those boundaries). 0 disables relaxation.");
  params.addParam<unsigned int>(
      "relaxation_iterations", 10, "Smoothing sweeps over the free nodes per relaxed move");
  params.addParam<std::vector<SubdomainName>>(
      "block", "Only elements in these subdomains can be moved or assigned to the crack sides");
  params.addRequiredParam<subdomain_id_type>(
      "lower_block_id", "Subdomain for crack elements on the side opposite to the normal");
  params.addRequiredParam<subdomain_id_type>(
      "upper_block_id", "Subdomain for crack elements on the side the normal points to");
  params.addParam<SubdomainName>("lower_block_name", "Name of the 'lower_block_id' subdomain");
  params.addParam<SubdomainName>("upper_block_name", "Name of the 'upper_block_id' subdomain");
  params.addClassDescription(
      "Fits a plane to the ridge of a phase-field crack, moves nodes onto it so element faces tile "
      "the crack, and assigns the elements on its two sides to two subdomains for "
      "BreakMeshByBlockGenerator to split.");
  return params;
}

PlanarCrackGenerator::PlanarCrackGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _threshold(getParam<Real>("threshold")),
    _ridge_search_distance(getParam<Real>("ridge_search_distance")),
    _min_scaled_jacobian(getParam<Real>("min_scaled_jacobian")),
    _min_boundary_alignment(getParam<Real>("min_boundary_alignment")),
    _relaxation_layers(getParam<unsigned int>("relaxation_layers")),
    _relaxation_iterations(getParam<unsigned int>("relaxation_iterations")),
    _lower_block_id(getParam<subdomain_id_type>("lower_block_id")),
    _upper_block_id(getParam<subdomain_id_type>("upper_block_id"))
{
  if (isParamValid("normal") && getParam<RealVectorValue>("normal").norm() == 0)
    paramError("normal", "The crack normal must be nonzero.");
  if (_lower_block_id == _upper_block_id)
    paramError("upper_block_id", "'lower_block_id' and 'upper_block_id' must differ.");

  declareMeshProperty<std::vector<Point>>("moved_node_positions", {});
  declareMeshProperty<std::vector<Point>>("moved_node_offsets", {});
}

Real
PlanarCrackGenerator::crackRidgeValue(PhaseFieldSolution & phase_field,
                                      const std::vector<Point> & points,
                                      const Real reach) const
{
  Real peak = -std::numeric_limits<Real>::max();
  for (const auto & p : points)
  {
    const Point on_plane = p - ((p - _plane_point) * _plane_normal) * _plane_normal;
    peak = std::max(peak, phase_field.findRidge(on_plane, _plane_normal, reach).second);
  }
  return peak;
}

std::unique_ptr<MeshBase>
PlanarCrackGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  if (!mesh->is_replicated())
    mooseError("PlanarCrackGenerator requires a replicated mesh.");
  if (!mesh->is_prepared())
    mesh->find_neighbors();

  for (const auto id : {_lower_block_id, _upper_block_id})
    if (MooseMeshUtils::hasSubdomainID(*mesh, id))
      paramError(id == _lower_block_id ? "lower_block_id" : "upper_block_id",
                 "Subdomain ",
                 id,
                 " already exists in the mesh; choose an unused id.");

  std::set<subdomain_id_type> restricted_blocks;
  if (isParamValid("block"))
  {
    const auto & names = getParam<std::vector<SubdomainName>>("block");
    const auto ids = MooseMeshUtils::getSubdomainIDs(*mesh, names);
    for (const auto i : index_range(ids))
      if (!MooseMeshUtils::hasSubdomainID(*mesh, ids[i]))
        paramError("block", "The subdomain '", names[i], "' was not found in the mesh.");
    restricted_blocks.insert(ids.begin(), ids.end());
  }
  const auto allowed = [&](const Elem & elem)
  { return restricted_blocks.empty() || restricted_blocks.count(elem.subdomain_id()); };

  for (const auto & elem : mesh->active_element_ptr_range())
    if (allowed(*elem) && elem->type() != TET4 && elem->type() != TET10)
      mooseError("PlanarCrackGenerator supports TET4 and TET10 elements only; element ",
                 elem->id(),
                 " is ",
                 Utility::enum_to_string(elem->type()),
                 ".");

  PhaseFieldSolution phase_field(*this);

  // 1. Band
  std::vector<const Elem *> band;
  std::vector<Point> centroids;
  for (const auto & elem : mesh->active_element_ptr_range())
  {
    if (!allowed(*elem))
      continue;
    const Point c = elem->vertex_average();
    if (phase_field.value(c) >= _threshold)
    {
      band.push_back(elem);
      centroids.push_back(c);
    }
  }
  if (band.empty())
  {
    mooseWarning("No element reached d >= ", _threshold, "; the mesh is returned unchanged.");
    return mesh;
  }
  Real band_size = 0;
  for (const auto * elem : band)
    band_size = std::max(band_size, elem->hmax());
  const Real reach = _ridge_search_distance > 0 ? _ridge_search_distance : 2 * band_size;

  // 2. Plane
  const bool fit_normal = !isParamValid("normal");
  RealVectorValue normal;
  if (fit_normal)
  {
    // The dominant eigenvector of the structure tensor is the normal up to sign
    RankTwoTensor structure;
    for (const auto i : index_range(band))
      structure +=
          band[i]->volume() * RankTwoTensor::selfOuterProduct(phase_field.gradient(centroids[i]));
    std::vector<Real> eigenvalues;
    RankTwoTensor eigenvectors;
    structure.symmetricEigenvaluesEigenvectors(eigenvalues, eigenvectors);
    normal = eigenvectors.column(2);
  }
  else
    normal = getParam<RealVectorValue>("normal").unit();

  // The plane passes through the weighted mean of the ridge points. Its normal comes from the
  // phase-field gradients rather than a best fit of the ridge points, which is ill-conditioned for
  // thin bodies: through a thin plate the ridge points spread about as much across the thickness as
  // the ridge wobbles, so the fit tilts the plane into the thickness.
  std::vector<Point> ridge(band.size());
  std::vector<Real> weight(band.size());
  Real total_weight = 0;
  Point plane_point;
  for (const auto i : index_range(band))
  {
    const auto [offset, peak] = phase_field.findRidge(centroids[i], normal, reach);
    ridge[i] = centroids[i] + offset * normal;
    weight[i] = peak;
    plane_point += peak * ridge[i];
    total_weight += peak;
  }
  plane_point /= total_weight;
  Real rms = 0;
  for (const auto i : index_range(band))
    rms += weight[i] * Utility::pow<2>((ridge[i] - plane_point) * normal);
  rms = std::sqrt(rms / total_weight);

  if (fit_normal)
  {
    // Make the arbitrary sign deterministic: the largest component points in its positive direction
    unsigned int largest = 0;
    for (unsigned int i = 1; i < LIBMESH_DIM; ++i)
      if (std::abs(normal(i)) > std::abs(normal(largest)))
        largest = i;
    if (normal(largest) < 0)
      normal = -normal;
  }
  _plane_point = plane_point;
  _plane_normal = normal;
  if (rms > band_size)
    mooseWarning("The ridge of the phase field is ",
                 rms,
                 " from the fitted plane on average (RMS), more than the largest band element size ",
                 band_size,
                 "; the crack is not planar and flattening it moves it away from the damage.");

  const Real tolerance = 1e-10 * band_size;
  const auto phi = [&](const Point & p) { return (p - _plane_point) * _plane_normal; };
  const auto inside = [&](const Elem & elem, const unsigned int side = libMesh::invalid_uint)
  {
    std::vector<Point> points;
    if (side == libMesh::invalid_uint)
    {
      points.push_back(elem.vertex_average());
      for (unsigned int v = 0; v < 4; ++v)
        points.push_back(elem.point(v));
    }
    else
    {
      Point centroid;
      for (const auto v : elem.nodes_on_side(side))
        if (v < 4)
        {
          points.push_back(elem.point(v));
          centroid += elem.point(v) / 3;
        }
      points.push_back(centroid);
    }
    return crackRidgeValue(phase_field, points, reach) >= _threshold;
  };

  std::unordered_map<dof_id_type, std::vector<const Elem *>> nodes_to_elems;
  MeshTools::build_nodes_to_elem_map(*mesh, nodes_to_elems);

  // 3.-4. Snap one end of every plane-crossing edge, first of the elements inside the crack extent,
  // then of crack elements that still cross the plane (at the edge of the extent)
  const auto straddles = [&](const Elem & elem)
  {
    Real lowest = std::numeric_limits<Real>::max(), highest = -lowest;
    for (unsigned int v = 0; v < 4; ++v)
    {
      lowest = std::min(lowest, phi(elem.point(v)));
      highest = std::max(highest, phi(elem.point(v)));
    }
    return lowest < -tolerance && highest > tolerance;
  };
  const auto crossing_edges = [&](const std::vector<const Elem *> & elems)
  {
    std::map<std::pair<dof_id_type, dof_id_type>, Real> crossing;
    for (const auto * elem : elems)
      for (const auto & edge : tet_edges)
      {
        const Real phi_a = phi(elem->point(edge[0])), phi_b = phi(elem->point(edge[1]));
        if (phi_a * phi_b < 0 && std::abs(phi_a) > tolerance && std::abs(phi_b) > tolerance)
        {
          auto a = elem->node_id(edge[0]), b = elem->node_id(edge[1]);
          crossing[{std::min(a, b), std::max(a, b)}] = std::min(std::abs(phi_a), std::abs(phi_b));
        }
      }
    std::vector<std::pair<dof_id_type, dof_id_type>> edges;
    for (const auto & [edge, distance] : crossing)
      edges.push_back(edge);
    std::stable_sort(edges.begin(),
                     edges.end(),
                     [&](const auto & e1, const auto & e2) { return crossing[e1] < crossing[e2]; });
    return edges;
  };

  std::vector<const Elem *> to_snap;
  for (const auto & elem : mesh->active_element_ptr_range())
    if (allowed(*elem) && straddles(*elem) && inside(*elem))
      to_snap.push_back(elem);

  std::map<dof_id_type, Point> original_positions;
  std::set<dof_id_type> moved_vertices;
  Real largest_move = 0;
  std::size_t n_refused_corner = 0, n_refused_boundary = 0, n_refused_inverted = 0,
              n_refused_quality = 0;
  // Vertices on external boundaries, which relaxation leaves in place
  std::set<dof_id_type> boundary_vertices;
  for (const auto & elem : mesh->active_element_ptr_range())
    for (const auto s : elem->side_index_range())
      if (!elem->neighbor_ptr(s))
        for (const auto v : elem->nodes_on_side(s))
          if (v < 4)
            boundary_vertices.insert(elem->node_id(v));

  // Unit normals of the distinct external boundary surfaces through a vertex. Faces within 20
  // degrees are one surface, e.g. the facets of a curved boundary, with their averaged normal.
  const auto boundary_normals = [&](const dof_id_type node_id)
  {
    std::vector<RealVectorValue> face_normals;
    if (!boundary_vertices.count(node_id))
      return face_normals;
    for (const auto * elem : nodes_to_elems[node_id])
      for (const auto s : elem->side_index_range())
      {
        if (elem->neighbor_ptr(s))
          continue;
        std::vector<Point> vertices;
        bool contains = false;
        for (const auto v : elem->nodes_on_side(s))
          if (v < 4)
          {
            vertices.push_back(elem->point(v));
            contains = contains || elem->node_id(v) == node_id;
          }
        if (!contains)
          continue;
        const auto face_normal = (vertices[1] - vertices[0]).cross(vertices[2] - vertices[0]).unit();
        bool known = false;
        for (auto & other : face_normals)
          if (face_normal * other.unit() > std::cos(20. * libMesh::pi / 180))
          {
            other += face_normal;
            known = true;
            break;
          }
        if (!known)
          face_normals.push_back(face_normal);
      }
    for (auto & face_normal : face_normals)
      face_normal = face_normal.unit();
    return face_normals;
  };

  // A displacement restricted to the boundaries through a vertex: within the surface for one
  // surface, along the line where two meet, none at a corner
  const auto within_boundaries = [&](const std::vector<RealVectorValue> & face_normals,
                                     const RealVectorValue & displacement)
  {
    if (face_normals.empty())
      return displacement;
    if (face_normals.size() == 1)
      return RealVectorValue(displacement - (displacement * face_normals[0]) * face_normals[0]);
    if (face_normals.size() == 2)
    {
      const auto line = face_normals[0].cross(face_normals[1]);
      return line.norm() > 0 ? RealVectorValue((displacement * line.unit()) * line.unit())
                             : RealVectorValue();
    }
    return RealVectorValue();
  };

  // The move bringing a vertex onto the plane, sliding within the boundaries it lies on
  const auto plane_move = [&](const dof_id_type node_id, const bool count) -> std::optional<Point>
  {
    const auto & node = mesh->node_ref(node_id);
    const auto face_normals = boundary_normals(node_id);
    RealVectorValue direction;
    if (face_normals.empty())
      direction = _plane_normal;
    else if (face_normals.size() == 1)
      direction = _plane_normal - (_plane_normal * face_normals[0]) * face_normals[0];
    else if (face_normals.size() == 2)
      direction = face_normals[0].cross(face_normals[1]);
    else
    {
      n_refused_corner += count;
      return {};
    }
    const Real alignment = direction.norm() > 0 ? direction.unit() * _plane_normal : 0;
    if (std::abs(alignment) < _min_boundary_alignment)
    {
      n_refused_boundary += count;
      return {};
    }
    return Point(node - (phi(node) / alignment) * direction.unit());
  };

  // Whether a region's elements are acceptable after a move: none inverted, and the worst no worse
  // than the quality limit or, if the region already was below it, than the worst before.
  const auto acceptable = [&](const std::map<const Elem *, Real> & before, const bool count)
  {
    Real worst_before = std::numeric_limits<Real>::max();
    Real worst_after = std::numeric_limits<Real>::max();
    for (const auto & [elem, quality] : before)
    {
      worst_before = std::min(worst_before, quality);
      worst_after = std::min(worst_after, scaledJacobian(*elem));
    }
    if (worst_after <= 0)
    {
      n_refused_inverted += count;
      return false;
    }
    if (worst_after < std::min(_min_scaled_jacobian, worst_before))
    {
      n_refused_quality += count;
      return false;
    }
    return true;
  };

  std::set<dof_id_type> relaxed_vertices;
  // Moves a vertex onto the plane. With relax, free vertices within a few element layers are
  // smoothed to make room; everything is reverted if an element still ends up unacceptable.
  const auto try_move = [&](const dof_id_type node_id, const bool relax, const bool count)
  {
    const auto target = plane_move(node_id, count);
    if (!target)
      return false;

    // Free vertices around the node: not on the plane, in the allowed blocks (boundary vertices only
    // move within their boundaries)
    std::vector<dof_id_type> free_vertices;
    if (relax)
    {
      std::set<dof_id_type> seen = {node_id}, frontier = {node_id};
      for (unsigned int layer = 0; layer < _relaxation_layers; ++layer)
      {
        std::set<dof_id_type> next;
        for (const auto n : frontier)
          for (const auto * elem : nodes_to_elems[n])
            if (allowed(*elem))
              for (unsigned int v = 0; v < 4; ++v)
              {
                const auto m = elem->node_id(v);
                if (seen.insert(m).second && std::abs(phi(mesh->point(m))) > tolerance)
                {
                  free_vertices.push_back(m);
                  next.insert(m);
                }
              }
        frontier = std::move(next);
      }
    }

    std::map<const Elem *, Real> before;
    std::map<dof_id_type, Point> old_positions = {{node_id, mesh->point(node_id)}};
    for (const auto * elem : nodes_to_elems[node_id])
      before.emplace(elem, scaledJacobian(*elem));
    for (const auto m : free_vertices)
    {
      old_positions.emplace(m, mesh->point(m));
      for (const auto * elem : nodes_to_elems[m])
        before.emplace(elem, scaledJacobian(*elem));
    }

    mesh->node_ref(node_id) = *target;
    if (relax)
      for (unsigned int iteration = 0; iteration < _relaxation_iterations; ++iteration)
        for (const auto m : free_vertices)
        {
          auto & free_node = mesh->node_ref(m);
          std::set<dof_id_type> neighbors;
          for (const auto * elem : nodes_to_elems[m])
            for (unsigned int v = 0; v < 4; ++v)
              if (elem->node_id(v) != m)
                neighbors.insert(elem->node_id(v));
          Point average;
          for (const auto n : neighbors)
            average += mesh->point(n) / Real(neighbors.size());
          // Boundary vertices smooth within their boundaries
          average = free_node + within_boundaries(boundary_normals(m), average - free_node);
          // Stay strictly on the node's side of the plane
          const Real side = phi(free_node);
          if (phi(average) * side <= 0 || std::abs(phi(average)) <= 0.1 * std::abs(side))
            continue;
          const auto elem_quality = [&]()
          {
            Real lowest = std::numeric_limits<Real>::max();
            for (const auto * elem : nodes_to_elems[m])
              lowest = std::min(lowest, scaledJacobian(*elem));
            return lowest;
          };
          const Real quality_before = elem_quality();
          const Point previous = free_node;
          free_node = average;
          if (elem_quality() < quality_before)
            free_node = previous;
        }

    if (!acceptable(before, count))
    {
      for (const auto & [n, position] : old_positions)
        mesh->node_ref(n) = position;
      return false;
    }

    largest_move = std::max(largest_move, (mesh->point(node_id) - old_positions[node_id]).norm());
    for (const auto & [n, position] : old_positions)
      if ((mesh->point(n) - position).norm() > tolerance)
      {
        original_positions.emplace(n, position);
        if (n == node_id)
        {
          moved_vertices.insert(n);
          relaxed_vertices.erase(n);
        }
        else if (!moved_vertices.count(n))
          relaxed_vertices.insert(n);
      }
    return true;
  };

  std::set<std::pair<dof_id_type, dof_id_type>> unresolved;
  struct PlaneFace
  {
    const Elem * elem;
    const Elem * neighbor;
    bool crack;
  };
  std::map<std::array<dof_id_type, 3>, PlaneFace> plane_faces;
  std::set<const Elem *> crack_elems;
  const unsigned int max_passes = 10;
  for (unsigned int pass = 0; pass < max_passes; ++pass)
  {
    const auto n_moved_before = moved_vertices.size();
    for (const auto & [a, b] : crossing_edges(to_snap))
    {
      const Real phi_a = phi(mesh->point(a)), phi_b = phi(mesh->point(b));
      // An earlier move may already have put one end on the plane
      if (std::abs(phi_a) <= tolerance || std::abs(phi_b) <= tolerance)
        continue;
      const auto nearer = std::abs(phi_a) <= std::abs(phi_b) ? a : b;
      const auto farther = nearer == a ? b : a;
      const bool relaxing = _relaxation_layers > 0;
      if (!try_move(nearer, false, !relaxing) && !try_move(farther, false, !relaxing) &&
          (!relaxing || (!try_move(nearer, true, true) && !try_move(farther, true, true))))
        unresolved.insert({a, b});
    }

    // 5. Crack faces: interior faces on the plane inside the crack extent
    std::set<dof_id_type> plane_vertices;
    for (const auto & elem : mesh->active_element_ptr_range())
      for (unsigned int v = 0; v < 4; ++v)
        if (std::abs(phi(elem->point(v))) <= tolerance)
          plane_vertices.insert(elem->node_id(v));

    // Straighten the edges of every element with a vertex on the plane, so faces on the plane are
    // planar also where their vertices were already on it and the edges were curved.
    std::set<const Elem *> straightened;
    for (const auto & vertices : {plane_vertices, relaxed_vertices})
      for (const auto node_id : vertices)
        for (const auto * elem : nodes_to_elems[node_id])
          if (elem->type() == TET10)
            straightened.insert(elem);
    for (const auto * const_elem : straightened)
    {
      auto * elem = mesh->elem_ptr(const_elem->id());
      for (const auto e : elem->edge_index_range())
      {
        const auto edge_nodes = elem->nodes_on_edge(e);
        auto & middle = elem->node_ref(edge_nodes[2]);
        const Point midpoint = 0.5 * (elem->point(edge_nodes[0]) + elem->point(edge_nodes[1]));
        if ((midpoint - middle).norm() > tolerance)
        {
          original_positions.emplace(middle.id(), middle);
          middle = midpoint;
        }
      }
    }

    plane_faces.clear();
    for (const auto node_id : plane_vertices)
      for (const auto * elem : nodes_to_elems[node_id])
        for (const auto s : elem->side_index_range())
        {
          const Elem * neighbor = elem->neighbor_ptr(s);
          if (!neighbor || neighbor == remote_elem)
            continue;
          std::array<dof_id_type, 3> key;
          unsigned int n_vertices = 0;
          bool on_plane = true;
          for (const auto v : elem->nodes_on_side(s))
            if (v < 4)
            {
              key[n_vertices++] = elem->node_id(v);
              on_plane = on_plane && plane_vertices.count(elem->node_id(v));
            }
          if (!on_plane)
            continue;
          std::sort(key.begin(), key.end());
          if (!plane_faces.count(key))
            plane_faces[key] = {elem, neighbor, inside(*elem, s)};
        }

    // Crack nodes are surrounded by crack faces; the others on crack faces form the front
    std::map<dof_id_type, bool> surrounded;
    for (const auto & [key, face] : plane_faces)
      for (const auto node_id : key)
      {
        auto it = surrounded.try_emplace(node_id, true).first;
        it->second = it->second && face.crack;
      }

    crack_elems.clear();
    for (const auto & [key, face] : plane_faces)
      if (face.crack)
      {
        crack_elems.insert(face.elem);
        crack_elems.insert(face.neighbor);
        for (const auto node_id : key)
          if (surrounded[node_id])
            for (const auto * elem : nodes_to_elems[node_id])
              crack_elems.insert(elem);
      }

    // Crack elements still crossing the plane get their edges snapped in the next pass
    to_snap.clear();
    for (const auto * elem : crack_elems)
      if (allowed(*elem) && straddles(*elem))
        to_snap.push_back(elem);
    if (to_snap.empty() || moved_vertices.size() == n_moved_before)
      break;
  }
  std::size_t n_unresolved = 0;
  for (const auto & [a, b] : unresolved)
    n_unresolved += phi(mesh->point(a)) * phi(mesh->point(b)) < 0 &&
                    std::abs(phi(mesh->point(a))) > tolerance &&
                    std::abs(phi(mesh->point(b))) > tolerance;

  std::set<const Elem *> touched;
  for (const auto node_id : moved_vertices)
    for (const auto * elem : nodes_to_elems[node_id])
      touched.insert(elem);

  // Final sides. Seeds: both elements of each crack face on the plane, and elements inside the
  // extent that still cross the plane because their edges could not be snapped. Every element goes
  // to the side its centroid is on: for snapped elements that is the side of all their vertices,
  // and around unsnapped ones the crack follows their faces (sawtooth) instead of leaving a hole.
  // Then, as in PhaseFieldCrackSubdomainGenerator, every node between the two sides where the ridge
  // reaches the threshold takes all elements around it, so BreakMeshByBlockGenerator duplicates it;
  // nodes where the ridge is below the threshold stay shared as the crack front.
  std::map<const Elem *, bool> is_upper;
  std::size_t n_outside_blocks = 0;
  const auto label = [&](const Elem * elem)
  {
    if (is_upper.count(elem))
      return false;
    if (!allowed(*elem))
    {
      ++n_outside_blocks;
      return false;
    }
    is_upper[elem] = phi(elem->vertex_average()) >= 0;
    return true;
  };
  for (const auto & [key, face] : plane_faces)
    if (face.crack)
    {
      label(face.elem);
      label(face.neighbor);
    }
  std::size_t n_sawtooth_elems = 0;
  for (const auto & elem : mesh->active_element_ptr_range())
    if (allowed(*elem) && straddles(*elem) && inside(*elem))
      n_sawtooth_elems += label(elem);

  std::set<dof_id_type> checked_nodes;
  bool added = true;
  while (added)
  {
    added = false;
    std::set<dof_id_type> interface_nodes;
    for (const auto & [elem, upper] : is_upper)
      for (const auto s : elem->side_index_range())
      {
        const Elem * neighbor = elem->neighbor_ptr(s);
        if (!neighbor || neighbor == remote_elem)
          continue;
        const auto it = is_upper.find(neighbor);
        if (it != is_upper.end() && it->second != upper)
          for (const auto v : elem->nodes_on_side(s))
            interface_nodes.insert(elem->node_id(v));
      }
    for (const auto node_id : interface_nodes)
    {
      if (!checked_nodes.insert(node_id).second)
        continue;
      // The ridge through the node and the centroids around it, which smooths out single dips
      std::vector<Point> points = {mesh->point(node_id)};
      for (const auto * elem : nodes_to_elems[node_id])
        points.push_back(elem->vertex_average());
      if (crackRidgeValue(phase_field, points, reach) < _threshold)
        continue;
      for (const auto * elem : nodes_to_elems[node_id])
        added = label(elem) || added;
    }
  }

  std::size_t n_lower = 0, n_upper = 0, n_crack_faces = 0, n_off_plane_faces = 0;
  for (const auto & [elem, upper] : is_upper)
  {
    mesh->elem_ptr(elem->id())->subdomain_id() = upper ? _upper_block_id : _lower_block_id;
    ++(upper ? n_upper : n_lower);
    if (!upper)
      for (const auto s : elem->side_index_range())
      {
        const auto it = is_upper.find(elem->neighbor_ptr(s));
        if (it == is_upper.end() || !it->second)
          continue;
        ++n_crack_faces;
        for (const auto v : elem->nodes_on_side(s))
          if (v < 4 && std::abs(phi(elem->point(v))) > tolerance)
          {
            ++n_off_plane_faces;
            break;
          }
      }
  }
  if (isParamValid("lower_block_name"))
    mesh->subdomain_name(_lower_block_id) = getParam<SubdomainName>("lower_block_name");
  if (isParamValid("upper_block_name"))
    mesh->subdomain_name(_upper_block_id) = getParam<SubdomainName>("upper_block_name");

  // 6. Moved nodes for looking fields up at their original positions
  auto & positions = setMeshProperty<std::vector<Point>>("moved_node_positions");
  auto & offsets = setMeshProperty<std::vector<Point>>("moved_node_offsets");
  positions.clear();
  offsets.clear();
  for (const auto & [node_id, original] : original_positions)
  {
    positions.push_back(mesh->point(node_id));
    offsets.push_back(original - mesh->point(node_id));
  }

  Real lowest_quality = std::numeric_limits<Real>::max();
  for (const auto * elem : touched)
    lowest_quality = std::min(lowest_quality, scaledJacobian(*elem));

  _console << name() << ": plane through " << _plane_point << " with normal " << _plane_normal
           << " (ridge RMS distance " << rms << ", " << band.size() << " band elements); moved "
           << moved_vertices.size() << " vertices (largest move " << largest_move << ") and "
           << relaxed_vertices.size() << " relaxed vertices and "
           << original_positions.size() - moved_vertices.size() - relaxed_vertices.size()
           << " mid-edge nodes"
           << (touched.empty() ? "" : "; lowest scaled Jacobian of moved elements ")
           << (touched.empty() ? "" : std::to_string(lowest_quality)) << "; " << n_crack_faces
           << " crack faces (" << n_off_plane_faces << " off the plane); " << n_lower
           << " lower and " << n_upper << " upper elements." << std::endl;
  if (n_unresolved)
    mooseWarning(n_unresolved,
                 " plane-crossing edges could not be snapped; around the ",
                 n_sawtooth_elems,
                 " elements still crossing the plane the crack follows element faces (",
                 n_off_plane_faces,
                 " crack faces off the plane). Refused moves: ",
                 n_refused_corner,
                 " boundary corners, ",
                 n_refused_boundary,
                 " boundary slides not aligned with the normal ('min_boundary_alignment'), ",
                 n_refused_inverted,
                 " inverting an element, ",
                 n_refused_quality,
                 " below 'min_scaled_jacobian'.");
  if (n_outside_blocks)
    mooseWarning(n_outside_blocks,
                 " elements needed on the crack sides are outside 'block' and were not assigned; "
                 "the crack stays closed there.");

  mesh->unset_is_prepared();
  return mesh;
}
