#include "PlanarCrackGenerator.h"
#include "MooseMeshUtils.h"
#include "PhaseFieldSolution.h"
#include "RankTwoTensor.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/remote_elem.h"
#include "libmesh/string_to_enum.h"

#include <numeric>
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
  params.addParam<bool>(
      "bisect_unsnapped_edges",
      true,
      "Split every edge that still crosses the plane after snapping at the crossing point (all "
      "elements around it are bisected, so the mesh stays conforming), making the crack exactly "
      "planar. Without it the crack follows element faces around those edges.");
  params.addRangeCheckedParam<Real>(
      "existing_crack_angle",
      20,
      "existing_crack_angle >= 0 & existing_crack_angle < 90",
      "External boundary faces inside the crack extent whose normal is within this angle (degrees) "
      "of the crack normal, and within 'existing_crack_opening' of the plane, are the faces of a "
      "crack cut earlier. Their nodes are not moved and the elements behind them are assigned to "
      "the crack sides, so a mesh that already contains a crack keeps it.");
  params.addRangeCheckedParam<Real>(
      "existing_crack_opening",
      0,
      "existing_crack_opening >= 0",
      "Largest distance from the plane at which a boundary face counts as an existing crack face, "
      "i.e. half the widest crack opening expected. 0 uses the largest band element size.");
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
    _bisect_unsnapped_edges(getParam<bool>("bisect_unsnapped_edges")),
    _existing_crack_angle(getParam<Real>("existing_crack_angle")),
    _existing_crack_opening(getParam<Real>("existing_crack_opening")),
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
    // Machine-readable summary, parsed by scripts deciding whether a crack was cut
    _console << name() << ": crack summary: band_elements=0 crack_faces=0 existing_faces=0"
             << std::endl;
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

  // Faces of a crack cut in an earlier step: external boundary faces in the crack plane facing
  // along the crack normal. Their nodes stay where they are and the elements behind them are
  // assigned to the crack sides below, so the existing crack is kept and handed on. The search
  // starts from faces inside the extent and spreads along the crack surface, so parts where the
  // damage has since dropped below the threshold are kept too.
  const Real opening = _existing_crack_opening > 0 ? _existing_crack_opening : band_size;
  const Real min_alignment = std::cos(_existing_crack_angle * libMesh::pi / 180);
  std::vector<const Elem *> existing_crack_elems;
  std::set<dof_id_type> frozen_vertices;
  const auto find_existing_crack = [&]()
  {
    existing_crack_elems.clear();
    frozen_vertices.clear();

    struct CrackFace
    {
      const Elem * elem;
      unsigned int side;
      std::vector<dof_id_type> vertices;
      bool seed;
    };
    std::vector<CrackFace> candidates;
    for (const auto & elem : mesh->active_element_ptr_range())
    {
      if (!allowed(*elem))
        continue;
      for (const auto s : elem->side_index_range())
      {
        if (elem->neighbor_ptr(s))
          continue;
        std::vector<Point> points;
        std::vector<dof_id_type> vertices;
        Real furthest = 0;
        for (const auto v : elem->nodes_on_side(s))
          if (v < 4)
          {
            points.push_back(elem->point(v));
            vertices.push_back(elem->node_id(v));
            furthest = std::max(furthest, std::abs(phi(elem->point(v))));
          }
        if (furthest > opening)
          continue;
        const auto face_normal = (points[1] - points[0]).cross(points[2] - points[0]).unit();
        if (std::abs(face_normal * _plane_normal) < min_alignment)
          continue;
        std::sort(vertices.begin(), vertices.end());
        candidates.push_back({elem, s, vertices, inside(*elem, s)});
      }
    }

    // Spread from the seeds over candidates sharing an edge, i.e. two vertices
    std::map<std::pair<dof_id_type, dof_id_type>, std::vector<std::size_t>> by_edge;
    for (const auto i : index_range(candidates))
      for (const auto & [a, b] : {std::make_pair(0, 1), std::make_pair(1, 2), std::make_pair(0, 2)})
        by_edge[{candidates[i].vertices[a], candidates[i].vertices[b]}].push_back(i);

    std::vector<bool> on_crack(candidates.size(), false);
    std::queue<std::size_t> queue;
    for (const auto i : index_range(candidates))
      if (candidates[i].seed)
      {
        on_crack[i] = true;
        queue.push(i);
      }
    while (!queue.empty())
    {
      const auto i = queue.front();
      queue.pop();
      for (const auto & [a, b] : {std::make_pair(0, 1), std::make_pair(1, 2), std::make_pair(0, 2)})
        for (const auto j : by_edge[{candidates[i].vertices[a], candidates[i].vertices[b]}])
          if (!on_crack[j])
          {
            on_crack[j] = true;
            queue.push(j);
          }
    }

    for (const auto i : index_range(candidates))
      if (on_crack[i])
      {
        existing_crack_elems.push_back(candidates[i].elem);
        for (const auto v : candidates[i].elem->nodes_on_side(candidates[i].side))
          frozen_vertices.insert(candidates[i].elem->node_id(v));
      }
  };

  find_existing_crack();
  if (!existing_crack_elems.empty())
  {
    // Anchor the plane to the existing crack so new growth continues it instead of stepping off it
    Point anchor;
    RealVectorValue anchor_normal;
    Real total_area = 0;
    for (const auto * elem : existing_crack_elems)
      for (const auto s : elem->side_index_range())
      {
        if (elem->neighbor_ptr(s))
          continue;
        std::vector<Point> points;
        for (const auto v : elem->nodes_on_side(s))
          if (v < 4)
            points.push_back(elem->point(v));
        if (points.size() != 3)
          continue;
        const auto cross = (points[1] - points[0]).cross(points[2] - points[0]);
        const Real area = 0.5 * cross.norm();
        if (area == 0 || std::abs(cross.unit() * _plane_normal) < min_alignment ||
            std::abs(phi(points[0])) > opening)
          continue;
        const auto face_normal = cross.unit() * _plane_normal < 0 ? -cross.unit() : cross.unit();
        anchor += area * (points[0] + points[1] + points[2]) / 3;
        anchor_normal += area * face_normal;
        total_area += area;
      }
    if (total_area > 0 && anchor_normal.norm() > 0)
    {
      const auto anchored = anchor_normal.unit();
      const Real turn = std::acos(std::min(1., std::abs(anchored * _plane_normal))) * 180 /
                        libMesh::pi;
      if (turn > _existing_crack_angle)
        mooseWarning("The phase field's ridge runs ",
                     turn,
                     " degrees away from the existing crack's plane, more than "
                     "'existing_crack_angle'. The crack is extended along the existing plane, so a "
                     "crack that turns is flattened onto it.");
      _plane_point = anchor / total_area;
      _plane_normal = anchored;
      rms = 0;
      for (const auto i : index_range(band))
        rms += weight[i] * Utility::pow<2>((ridge[i] - _plane_point) * _plane_normal);
      rms = std::sqrt(rms / total_weight);
      // The plane moved, so the faces on it and the extent are found again
      find_existing_crack();
    }
  }

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
              n_refused_quality = 0, n_refused_frozen = 0;
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
    // Nodes of an existing crack's faces keep the geometry that crack was cut with
    if (frozen_vertices.count(node_id))
    {
      n_refused_frozen += count;
      return false;
    }
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
                if (seen.insert(m).second && !frozen_vertices.count(m) &&
                    std::abs(phi(mesh->point(m))) > tolerance)
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
  // Splits every edge of the given elements that crosses the plane at the crossing point. Each
  // element around such an edge is replaced by two, one with each edge end swapped for the new node,
  // so the mesh stays conforming. The new node lies on the plane, so no edge through it crosses it,
  // and every split removes one crossing. Returns the number of edges split.
  const auto edge_key = [](const dof_id_type a, const dof_id_type b)
  { return std::make_pair(std::min(a, b), std::max(a, b)); };
  const auto original_of = [&](const dof_id_type node_id) -> Point
  {
    const auto it = original_positions.find(node_id);
    return it == original_positions.end() ? Point(mesh->point(node_id)) : it->second;
  };
  auto & boundary_info = mesh->get_boundary_info();
  std::size_t n_bisected = 0, n_split_elems = 0;
  const auto bisect = [&](const std::vector<const Elem *> & elems) -> std::size_t
  {
    // Keyed and ordered by id, never by address, so every processor of a replicated mesh splits in
    // the same order and gives new nodes and elements the same ids
    std::map<dof_id_type, std::set<dof_id_type>> vertex_elems;
    std::map<std::pair<dof_id_type, dof_id_type>, dof_id_type> edge_mid;
    for (auto & elem : mesh->active_element_ptr_range())
    {
      for (unsigned int v = 0; v < 4; ++v)
        vertex_elems[elem->node_id(v)].insert(elem->id());
      if (elem->type() == TET10)
        for (const auto e : elem->edge_index_range())
        {
          const auto nodes = elem->nodes_on_edge(e);
          edge_mid[edge_key(elem->node_id(nodes[0]), elem->node_id(nodes[1]))] =
              elem->node_id(nodes[2]);
        }
    }
    // The mid-edge node of an edge, created at its midpoint if the edge is new. New nodes take the
    // processor of the element they are made for; an unset one makes libMesh wait for its owner.
    processor_id_type owner = 0;
    const auto mid_node = [&](const dof_id_type a, const dof_id_type b) -> Node *
    {
      const auto key = edge_key(a, b);
      const auto it = edge_mid.find(key);
      if (it != edge_mid.end())
        return mesh->node_ptr(it->second);
      Node * mid =
          mesh->add_point(0.5 * (mesh->point(a) + mesh->point(b)), DofObject::invalid_id, owner);
      const Point original = 0.5 * (original_of(a) + original_of(b));
      if ((original - *mid).norm() > tolerance)
        original_positions.emplace(mid->id(), original);
      edge_mid[key] = mid->id();
      return mid;
    };
    const auto set_vertices = [&](Elem * elem, const std::array<dof_id_type, 4> & vertices)
    {
      for (unsigned int k = 0; k < 4; ++k)
        elem->set_node(k, mesh->node_ptr(vertices[k]));
      if (elem->type() == TET10)
        for (const auto e : elem->edge_index_range())
        {
          const auto nodes = elem->nodes_on_edge(e);
          elem->set_node(nodes[2], mid_node(vertices[nodes[0]], vertices[nodes[1]]));
        }
    };

    std::set<std::pair<dof_id_type, dof_id_type>> edges;
    for (const auto & [a, b] : crossing_edges(elems))
      edges.insert(edge_key(a, b));

    std::size_t n_split = 0;
    std::set<dof_id_type> replaced_mids;
    for (const auto & [a, b] : edges)
    {
      const Real phi_a = phi(mesh->point(a)), phi_b = phi(mesh->point(b));
      if (!(phi_a * phi_b < 0 && std::abs(phi_a) > tolerance && std::abs(phi_b) > tolerance))
        continue;
      std::vector<Elem *> around;
      for (const auto elem_id : vertex_elems[a])
        if (vertex_elems[b].count(elem_id))
          around.push_back(mesh->elem_ptr(elem_id));
      if (around.empty())
        continue;

      const Real t = phi_a / (phi_a - phi_b);
      Point crossing = mesh->point(a) + t * (mesh->point(b) - mesh->point(a));
      crossing -= phi(crossing) * _plane_normal;
      owner = around.front()->processor_id();
      Node * m = mesh->add_point(crossing, DofObject::invalid_id, owner);
      const Point original = original_of(a) + t * (original_of(b) - original_of(a));
      if ((original - crossing).norm() > tolerance)
        original_positions.emplace(m->id(), original);
      const auto old_mid = edge_mid.find(edge_key(a, b));
      if (old_mid != edge_mid.end())
      {
        replaced_mids.insert(old_mid->second);
        edge_mid.erase(old_mid);
      }

      for (auto * parent : around)
      {
        owner = parent->processor_id();
        const unsigned int ia = parent->local_node(a), ib = parent->local_node(b);
        std::array<dof_id_type, 4> first, second;
        for (unsigned int k = 0; k < 4; ++k)
          first[k] = second[k] = parent->node_id(k);
        first[ib] = m->id();
        second[ia] = m->id();

        // The side opposite a vertex is the one that does not contain it
        const auto opposite = [&](const unsigned int k) -> unsigned int
        {
          for (const auto side : parent->side_index_range())
            if (!parent->is_node_on_side(k, side))
              return side;
          return libMesh::invalid_uint;
        };
        const unsigned int side_opposite_a = opposite(ia), side_opposite_b = opposite(ib);
        std::vector<std::vector<boundary_id_type>> side_ids(parent->n_sides());
        for (const auto side : parent->side_index_range())
          boundary_info.boundary_ids(parent, side, side_ids[side]);

        auto child = Elem::build(parent->type());
        child->subdomain_id() = parent->subdomain_id();
        child->processor_id() = parent->processor_id();
        set_vertices(child.get(), second);
        Elem * added = mesh->add_elem(std::move(child));
        for (const auto i : make_range(parent->n_extra_integers()))
          added->set_extra_integer(i, parent->get_extra_integer(i));
        set_vertices(parent, first);

        // Each half keeps the parts of the parent's faces it contains. The face through the new node
        // and the edge opposite the split one is new and interior.
        boundary_info.remove_side(parent, side_opposite_a);
        for (const auto side : parent->side_index_range())
          if (side != side_opposite_b)
            for (const auto id : side_ids[side])
              boundary_info.add_side(added, side, id);

        vertex_elems[b].erase(parent->id());
        vertex_elems[m->id()].insert(parent->id());
        for (unsigned int k = 0; k < 4; ++k)
          vertex_elems[second[k]].insert(added->id());
        ++n_split_elems;
      }
      ++n_split;
    }

    if (n_split)
    {
      mesh->find_neighbors();
      nodes_to_elems.clear();
      MeshTools::build_nodes_to_elem_map(*mesh, nodes_to_elems);
      // Mid-edge nodes of split edges are no longer used. They are left in place, since deleting
      // nodes leaves gaps in the node ids that BreakMeshByBlockGenerator mishandles; mesh
      // preparation removes unused nodes.
      for (const auto id : replaced_mids)
        if (!nodes_to_elems.count(id))
          original_positions.erase(id);
    }
    return n_split;
  };

  // Snap and split, assign the sides, and repeat while the side assignment took elements that
  // still cross the plane (e.g. around crack nodes at the edge of the extent), so every crack face
  // ends up on the plane.
  std::map<const Elem *, bool> is_upper;
  std::size_t n_outside_blocks = 0, n_sawtooth_elems = 0, n_pinched_nodes = 0;
  const unsigned int max_rounds = 5;
  for (unsigned int round = 0; round < max_rounds; ++round)
  {
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

      // Split what snapping left crossing the plane, so the crack is exactly planar
      std::size_t n_bisected_now = 0;
      if (_bisect_unsnapped_edges)
      {
        std::vector<const Elem *> crossing_elems;
        for (const auto * elem : to_snap)
          if (straddles(*elem))
            crossing_elems.push_back(elem);
        n_bisected_now = bisect(crossing_elems);
        n_bisected += n_bisected_now;
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
      if (to_snap.empty() || (moved_vertices.size() == n_moved_before && n_bisected_now == 0))
        break;
    }

    // Final sides. Seeds: both elements of each crack face on the plane, and elements inside the
    // extent that still cross the plane because their edges could not be snapped. Every element goes
    // to the side its centroid is on: for snapped elements that is the side of all their vertices,
    // and around unsnapped ones the crack follows their faces (sawtooth) instead of leaving a hole.
    // Then, as in PhaseFieldCrackSubdomainGenerator, every node between the two sides where the ridge
    // reaches the threshold takes all elements around it, so BreakMeshByBlockGenerator duplicates it;
    // nodes where the ridge is below the threshold stay shared as the crack front.
    is_upper.clear();
    n_outside_blocks = 0;
    n_sawtooth_elems = 0;
    n_pinched_nodes = 0;
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
    for (const auto & elem : mesh->active_element_ptr_range())
      if (allowed(*elem) && straddles(*elem) && inside(*elem))
        n_sawtooth_elems += label(elem);
    // The elements behind an existing crack's faces, so its two sides are subdomains again
    for (const auto * elem : existing_crack_elems)
      label(elem);

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
      if (added)
        continue;

      // The crack surface must be a manifold for BreakMeshByBlockGenerator's output to be a valid
      // solid. A node that stays shared must have its crack faces in one fan (faces joined across
      // edges to open nodes) that ends at a front edge somewhere (an edge of only one crack face,
      // not on the external boundary). Otherwise the crack is pinched there, e.g. along an edge
      // with crack faces on both sides whose nodes stayed shared because the ridge dips below the
      // threshold, and the node is opened as well.
      std::map<dof_id_type, std::vector<std::array<dof_id_type, 3>>> crack_faces_of;
      for (const auto & [elem, upper] : is_upper)
        for (const auto s : elem->side_index_range())
        {
          const Elem * neighbor = elem->neighbor_ptr(s);
          if (!neighbor || neighbor == remote_elem || neighbor->id() < elem->id())
            continue;
          const auto it = is_upper.find(neighbor);
          if (it == is_upper.end() || it->second == upper)
            continue;
          std::array<dof_id_type, 3> key;
          unsigned int n_vertices = 0;
          for (const auto v : elem->nodes_on_side(s))
            if (v < 4)
              key[n_vertices++] = elem->node_id(v);
          for (const auto node_id : key)
            crack_faces_of[node_id].push_back(key);
        }
      const auto open = [&](dof_id_type node_id)
      {
        for (const auto * elem : nodes_to_elems[node_id])
          if (!is_upper.count(elem))
            return false;
        return true;
      };
      const auto on_external_boundary = [&](dof_id_type a, dof_id_type b)
      {
        for (const auto * elem : nodes_to_elems[a])
          for (const auto s : elem->side_index_range())
            if (!elem->neighbor_ptr(s))
            {
              unsigned int n_found = 0;
              for (const auto v : elem->nodes_on_side(s))
                n_found += v < 4 && (elem->node_id(v) == a || elem->node_id(v) == b);
              if (n_found == 2)
                return true;
            }
        return false;
      };
      std::vector<dof_id_type> pinched;
      for (const auto & [node_id, faces] : crack_faces_of)
      {
        if (open(node_id))
          continue;
        std::map<dof_id_type, std::vector<std::size_t>> edge_faces;
        for (const auto f : index_range(faces))
          for (const auto other : faces[f])
            if (other != node_id)
              edge_faces[other].push_back(f);
        std::vector<std::size_t> parent(faces.size());
        std::iota(parent.begin(), parent.end(), 0);
        const auto find = [&](std::size_t f)
        {
          while (parent[f] != f)
            f = parent[f] = parent[parent[f]];
          return f;
        };
        bool front_end = false;
        for (const auto & [other, edge] : edge_faces)
          if (edge.size() == 1)
            front_end = front_end || !on_external_boundary(node_id, other);
          else if (open(other))
            for (const auto f : edge)
              parent[find(f)] = find(edge.front());
        std::set<std::size_t> fans;
        for (const auto f : index_range(faces))
          fans.insert(find(f));
        if (fans.size() > 1 || !front_end)
          pinched.push_back(node_id);
      }
      for (const auto node_id : pinched)
        for (const auto * elem : nodes_to_elems[node_id])
          added = label(elem) || added;
      n_pinched_nodes += pinched.size();
    }

    // Elements assigned to a side that still cross the plane are snapped or split in another round
    to_snap.clear();
    for (const auto & [elem, upper] : is_upper)
      if (straddles(*elem))
        to_snap.push_back(elem);
    if (to_snap.empty() || !_bisect_unsnapped_edges)
      break;
  }

  // Edges that still exist (not split) and still cross the plane
  std::size_t n_unresolved = 0;
  for (const auto & [a, b] : unresolved)
  {
    bool exists = false;
    for (const auto * elem : nodes_to_elems[a])
      for (unsigned int v = 0; v < 4; ++v)
        exists = exists || elem->node_id(v) == b;
    n_unresolved += exists && phi(mesh->point(a)) * phi(mesh->point(b)) < 0 &&
                    std::abs(phi(mesh->point(a))) > tolerance &&
                    std::abs(phi(mesh->point(b))) > tolerance;
  }

  std::set<const Elem *> touched;
  for (const auto node_id : moved_vertices)
    for (const auto * elem : nodes_to_elems[node_id])
      touched.insert(elem);

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
           << (n_bisected ? "; split " + std::to_string(n_bisected) + " edges crossing the plane (" +
                                std::to_string(n_split_elems) + " elements bisected)"
                          : "")
           << (touched.empty() ? "" : "; lowest scaled Jacobian of moved elements ")
           << (touched.empty() ? "" : std::to_string(lowest_quality)) << "; " << n_crack_faces
           << " crack faces (" << n_off_plane_faces << " off the plane); " << n_lower
           << " lower and " << n_upper << " upper elements"
           << (n_pinched_nodes ? "; opened " + std::to_string(n_pinched_nodes) +
                                     " crack front nodes that pinched the crack surface"
                               : "")
           << (existing_crack_elems.empty()
                   ? "."
                   : "; kept " + std::to_string(existing_crack_elems.size()) +
                         " faces of an existing crack.")
           << std::endl;
  // Machine-readable summary, parsed by scripts deciding whether a crack was cut
  _console << name() << ": crack summary: band_elements=" << band.size()
           << " crack_faces=" << n_crack_faces
           << " existing_faces=" << existing_crack_elems.size() << std::endl;
  if (n_crack_faces == 0 && existing_crack_elems.empty())
    mooseWarning("No crack faces were found where the ridge reaches d >= ",
                 _threshold,
                 "; no elements were assigned to the crack subdomains.");
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
                 " below 'min_scaled_jacobian', ",
                 n_refused_frozen,
                 " nodes of an existing crack.");
  if (n_outside_blocks)
    mooseWarning(n_outside_blocks,
                 " elements needed on the crack sides are outside 'block' and were not assigned; "
                 "the crack stays closed there.");

  mesh->unset_is_prepared();
  return mesh;
}
