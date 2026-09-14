#include "PhaseFieldCrackSubdomainGenerator.h"
#include "Conversion.h"
#include "KDTree.h"
#include "MooseMeshUtils.h"
#include "RankTwoTensor.h"

#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/remote_elem.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/string_to_enum.h"

#include <queue>

registerMooseObject("raccoonApp", PhaseFieldCrackSubdomainGenerator);

InputParameters
PhaseFieldCrackSubdomainGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");

  params.addRequiredParam<FileName>("solution_file",
                                    "The Exodus file holding the phase-field variable. Its mesh "
                                    "may differ from the input mesh; values are interpolated.");
  params.addParam<std::string>("variable", "d", "The nodal phase-field variable in the file");
  params.addParam<std::string>(
      "timestep", "LATEST", "The Exodus time step (1-based) to read, or 'LATEST'");
  params.addParam<MooseEnum>(
      "variable_order",
      MooseEnum("AUTO FIRST SECOND", "AUTO"),
      "Interpolation order of the phase-field variable. AUTO uses the order of the solution "
      "file's elements, so second-order meshes use their mid-edge values.");

  params.addRequiredRangeCheckedParam<Real>(
      "threshold",
      "threshold > 0 & threshold <= 1",
      "Elements whose phase-field value at the centroid is at least this form the crack band");

  params.addParam<RealVectorValue>(
      "normal",
      "A fixed crack normal. The side it points to becomes 'upper_block_id'. Use for cracks that "
      "are close to flat; give 'normal_radius' instead to estimate the normal per element.");
  params.addRangeCheckedParam<Real>(
      "normal_radius",
      "normal_radius > 0",
      "Radius of the neighborhood of band elements used to estimate the local crack normal. A "
      "value of about the band width works well. On each crack the side that becomes "
      "'upper_block_id' is arbitrary.");
  params.addRangeCheckedParam<Real>(
      "gradient_tolerance",
      0.05,
      "gradient_tolerance >= 0 & gradient_tolerance < 1",
      "Band elements with |grad(d) . n| below this fraction of the band's largest |grad(d)| are "
      "assigned from their neighbors instead of from their own gradient. Raise it if elements "
      "deep inside a d = 1 plateau land on the wrong side: interpolating the plateau edge leaves "
      "small spurious gradients there.");
  params.addParam<unsigned int>(
      "smoothing_passes",
      2,
      "Number of majority-vote passes; an element switches side when most of its band face "
      "neighbors are on the other side");

  params.addParam<std::vector<SubdomainName>>(
      "block", "Only elements in these subdomains can become part of the crack band");
  params.addRequiredParam<subdomain_id_type>(
      "lower_block_id", "Subdomain for band elements on the side opposite to the normal");
  params.addRequiredParam<subdomain_id_type>(
      "upper_block_id", "Subdomain for band elements on the side the normal points to");
  params.addParam<SubdomainName>("lower_block_name", "Name of the 'lower_block_id' subdomain");
  params.addParam<SubdomainName>("upper_block_name", "Name of the 'upper_block_id' subdomain");

  params.addClassDescription(
      "Assigns the elements of a phase-field crack band (d >= threshold) to two subdomains, one "
      "on each side of the crack mid-surface, for BreakMeshByBlockGenerator to split along.");
  return params;
}

PhaseFieldCrackSubdomainGenerator::PhaseFieldCrackSubdomainGenerator(
    const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _solution_file(getParam<FileName>("solution_file")),
    _variable(getParam<std::string>("variable")),
    _threshold(getParam<Real>("threshold")),
    _has_normal(isParamValid("normal")),
    _normal_radius(isParamValid("normal_radius") ? getParam<Real>("normal_radius") : 0),
    _gradient_tolerance(getParam<Real>("gradient_tolerance")),
    _smoothing_passes(getParam<unsigned int>("smoothing_passes")),
    _lower_block_id(getParam<subdomain_id_type>("lower_block_id")),
    _upper_block_id(getParam<subdomain_id_type>("upper_block_id"))
{
  if (_has_normal == isParamValid("normal_radius"))
    paramError("normal", "Exactly one of 'normal' and 'normal_radius' must be given.");
  if (_has_normal && getParam<RealVectorValue>("normal").norm() == 0)
    paramError("normal", "The crack normal must be nonzero.");
  if (_lower_block_id == _upper_block_id)
    paramError("upper_block_id", "'lower_block_id' and 'upper_block_id' must differ.");
}

std::unique_ptr<MeshBase>
PhaseFieldCrackSubdomainGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  // The side assignment walks the band across element neighbors, so every processor needs the
  // whole mesh to reach the same answer.
  if (!mesh->is_replicated())
    mooseError(name(), " requires a replicated mesh.");
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

  // Load the phase-field variable into a MeshFunction, as SolutionUserObject does.
  libMesh::ReplicatedMesh solution_mesh(comm());
  libMesh::ExodusII_IO exodus(solution_mesh);
  exodus.read(_solution_file);
  solution_mesh.allow_renumbering(false);
  solution_mesh.prepare_for_use();

  const auto & nodal_names = exodus.get_nodal_var_names();
  if (std::find(nodal_names.begin(), nodal_names.end(), _variable) == nodal_names.end())
    paramError("variable",
               "Nodal variable '",
               _variable,
               "' was not found in '",
               _solution_file,
               "'. Available: ",
               Moose::stringify(nodal_names));

  const int n_steps = exodus.get_num_time_steps();
  if (n_steps == 0)
    paramError("solution_file", "The file contains no time steps.");
  int timestep = n_steps;
  const auto & timestep_string = getParam<std::string>("timestep");
  if (timestep_string != "LATEST")
  {
    std::istringstream ss(timestep_string);
    if (!((ss >> timestep) && ss.eof()) || timestep < 1 || timestep > n_steps)
      paramError("timestep", "Expected 'LATEST' or an integer from 1 to ", n_steps, ".");
  }

  Order order = FIRST;
  const auto & order_enum = getParam<MooseEnum>("variable_order");
  if (order_enum == "AUTO")
  {
    for (const auto & elem : solution_mesh.active_element_ptr_range())
      order = std::max(order, elem->default_order());
  }
  else
    order = Utility::string_to_enum<Order>(order_enum);

  libMesh::EquationSystems es(solution_mesh);
  auto & system = es.add_system<libMesh::ExplicitSystem>("phase_field");
  const auto var_num = system.add_variable(_variable, order, LAGRANGE);
  es.init();
  exodus.copy_nodal_solution(system, _variable, _variable, timestep);

  auto serialized = NumericVector<Number>::build(comm());
  serialized->init(system.n_dofs(), false, SERIAL);
  system.solution->localize(*serialized);

  libMesh::MeshFunction phase_field(es, *serialized, system.get_dof_map(), var_num);
  phase_field.init();
  // Centroids outside the solution mesh read as undamaged.
  phase_field.enable_out_of_mesh_mode(Number(0));

  // Collect the band: elements whose centroid value reaches the threshold.
  std::vector<Elem *> band;
  std::vector<Point> centroids;
  std::vector<RealVectorValue> gradients;
  for (auto & elem : mesh->active_element_ptr_range())
  {
    if (!restricted_blocks.empty() && !restricted_blocks.count(elem->subdomain_id()))
      continue;
    const Point c = elem->vertex_average();
    if (phase_field(c) < _threshold)
      continue;
    band.push_back(elem);
    centroids.push_back(c);
    gradients.push_back(phase_field.gradient(c));
  }

  const auto n_band = band.size();
  if (n_band == 0)
  {
    mooseWarning(name(),
                 ": no element reached d >= ",
                 _threshold,
                 "; the mesh is returned unchanged.");
    return mesh;
  }

  // Face neighbors within the band
  std::unordered_map<dof_id_type, std::size_t> band_index;
  for (const auto i : make_range(n_band))
    band_index[band[i]->id()] = i;
  std::vector<std::vector<std::size_t>> band_neighbors(n_band);
  for (const auto i : make_range(n_band))
    for (const auto s : band[i]->side_index_range())
    {
      const Elem * neighbor = band[i]->neighbor_ptr(s);
      if (!neighbor || neighbor == remote_elem)
        continue;
      const auto it = band_index.find(neighbor->id());
      if (it != band_index.end())
        band_neighbors[i].push_back(it->second);
    }

  // Crack normal per band element; a zero vector marks an element whose normal is unknown.
  std::vector<RealVectorValue> normals(n_band);
  if (_has_normal)
  {
    const auto normal = getParam<RealVectorValue>("normal").unit();
    std::fill(normals.begin(), normals.end(), normal);
  }
  else
  {
    // The dominant eigenvector of sum(grad(d) grad(d)^T) over the neighborhood is the normal up to
    // sign. Being quadratic in grad(d), it is unaffected by grad(d) flipping across the ridge.
    KDTree tree(centroids, 10);
    std::vector<nanoflann::ResultItem<std::size_t, Real>> hits;
    std::vector<Real> largest_eigenvalue(n_band, 0);
    for (const auto i : make_range(n_band))
    {
      tree.radiusSearch(centroids[i], _normal_radius, hits);
      RankTwoTensor structure;
      for (const auto & hit : hits)
        structure += band[hit.first]->volume() * RankTwoTensor::selfOuterProduct(gradients[hit.first]);

      std::vector<Real> eigenvalues;
      RankTwoTensor eigenvectors;
      structure.symmetricEigenvaluesEigenvectors(eigenvalues, eigenvectors);
      largest_eigenvalue[i] = eigenvalues[2];
      normals[i] = eigenvectors.column(2);
    }

    // Neighborhoods with (almost) no gradient, e.g. inside a d = 1 plateau, have no usable normal.
    const Real max_eigenvalue =
        *std::max_element(largest_eigenvalue.begin(), largest_eigenvalue.end());
    for (const auto i : make_range(n_band))
      if (largest_eigenvalue[i] <= 1e-12 * max_eigenvalue)
        normals[i] = RealVectorValue();

    // Orient the normals consistently by walking the band, flipping each one to agree with the
    // element it was reached from. Elements without a normal inherit it. Each connected band
    // starts from an element that has a normal, so its overall sign is arbitrary.
    std::vector<bool> visited(n_band, false);
    for (const auto seed : make_range(n_band))
    {
      if (visited[seed] || normals[seed].norm_sq() == 0)
        continue;
      std::queue<std::size_t> queue;
      queue.push(seed);
      visited[seed] = true;
      while (!queue.empty())
      {
        const auto i = queue.front();
        queue.pop();
        for (const auto j : band_neighbors[i])
        {
          if (visited[j])
            continue;
          if (normals[j].norm_sq() == 0)
            normals[j] = normals[i];
          else if (normals[j] * normals[i] < 0)
            normals[j] = -normals[j];
          visited[j] = true;
          queue.push(j);
        }
      }
    }
  }

  // Side of the ridge: d increases along n below the ridge (lower) and decreases above it.
  const short lower = -1, upper = 1, undecided = 0;
  Real max_gradient = 0;
  for (const auto & g : gradients)
    max_gradient = std::max(max_gradient, g.norm());
  const Real tolerance = _gradient_tolerance * max_gradient;

  std::vector<short> side(n_band, undecided);
  for (const auto i : make_range(n_band))
  {
    const Real slope = gradients[i] * normals[i];
    if (std::abs(slope) > tolerance)
      side[i] = slope > 0 ? lower : upper;
  }

  // Fill undecided elements from decided face neighbors. Updating all elements of a sweep at
  // once advances both flanks at the same rate, so across a plateau they meet in the middle.
  std::size_t n_filled = 0;
  bool break_ties = false;
  while (true)
  {
    auto updated = side;
    bool changed = false;
    for (const auto i : make_range(n_band))
    {
      if (side[i] != undecided)
        continue;
      int vote = 0;
      short first_decided = undecided;
      for (const auto j : band_neighbors[i])
      {
        vote += side[j];
        if (first_decided == undecided)
          first_decided = side[j];
      }
      if (vote != 0)
        updated[i] = vote > 0 ? upper : lower;
      else if (break_ties && first_decided != undecided)
        updated[i] = first_decided;
      if (updated[i] != undecided)
      {
        changed = true;
        ++n_filled;
      }
    }
    side = std::move(updated);

    if (changed)
      break_ties = false;
    else if (!break_ties)
      break_ties = true;
    else
      break;
  }

  // Remaining undecided elements have no decided element anywhere in their part of the band.
  std::size_t n_isolated = 0;
  for (auto & s : side)
    if (s == undecided)
    {
      s = upper;
      ++n_isolated;
    }
  if (n_isolated)
    mooseWarning(name(),
                 ": ",
                 n_isolated,
                 " band elements had no usable phase-field gradient in their part of the band and "
                 "were assigned to the upper side.");

  // Majority-vote smoothing: switch sides when most band face neighbors are on the other side.
  std::size_t n_switched = 0;
  for (unsigned int pass = 0; pass < _smoothing_passes; ++pass)
  {
    auto updated = side;
    for (const auto i : make_range(n_band))
    {
      std::size_t disagree = 0;
      for (const auto j : band_neighbors[i])
        disagree += side[j] != side[i];
      if (2 * disagree > band_neighbors[i].size())
      {
        updated[i] = -side[i];
        ++n_switched;
      }
    }
    side = std::move(updated);
  }

  std::size_t n_lower = 0;
  for (const auto i : make_range(n_band))
  {
    band[i]->subdomain_id() = side[i] == lower ? _lower_block_id : _upper_block_id;
    n_lower += side[i] == lower;
  }
  if (isParamValid("lower_block_name"))
    mesh->subdomain_name(_lower_block_id) = getParam<SubdomainName>("lower_block_name");
  if (isParamValid("upper_block_name"))
    mesh->subdomain_name(_upper_block_id) = getParam<SubdomainName>("upper_block_name");

  _console << name() << ": " << n_band << " band elements (d >= " << _threshold << "): " << n_lower
           << " lower, " << n_band - n_lower << " upper; " << n_filled
           << " assigned from neighbors, " << n_switched << " switched by smoothing." << std::endl;

  mesh->unset_is_prepared();
  return mesh;
}
