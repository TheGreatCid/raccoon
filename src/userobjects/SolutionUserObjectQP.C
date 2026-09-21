#include "SolutionUserObjectQP.h"
#include "libmesh/nemesis_io.h"
// MOOSE includes
#include "ConsoleUtils.h"
#include "MooseError.h"
#include "MooseMesh.h"
#include "MooseUtils.h"
#include "MooseVariableFE.h"
#include "RotationMatrix.h"
#include "Function.h"
#include "MooseMeshUtils.h"

// libMesh includes
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/exodusII_io_helper.h"
#include "libmesh/enum_xdr_mode.h"
#include "Qp_Mapping.h"

#if defined(__GLIBC__)
#include <malloc.h>
#endif
registerMooseObject("raccoonApp", SolutionUserObjectQP);

InputParameters
SolutionUserObjectQP::validParams()
{
  InputParameters params = SolutionUserObject::validParams();
  params.addParam<std::vector<MaterialName>>("tensor_materials", "materials to output qps on");
  params.addParam<std::vector<MaterialName>>("materials", "materials to output qps on");
  params.addRequiredParam<MooseEnum>(
      "element", MooseEnum(QpMapping::ELEMENT_ENUM_DEFINITION), "The element type");
  params.addParam<bool>(
      "free_after_initial",
      false,
      "Release the loaded exodus solution/mesh payload (serialized-per-rank solution, "
      "mesh functions, equation systems, source mesh, and exodus reader) after the "
      "one-time initial recovery reads. Only safe when this object is queried "
      "exclusively during stateful-property initialization (the recovery pattern) and "
      "is not time-interpolating. Do not enable when the object is queried after "
      "INITIAL or when mesh adaptivity may re-initialize stateful properties.");
  params.addRangeCheckedParam<Real>(
      "relax_threshold",
      "relax_threshold > 0 & relax_threshold <= 1",
      "Relax the recovered state of damaged elements behind a cut crack: elements on the crack "
      "faces and those connected to them whose centroid value of 'relax_variable' reaches this, "
      "and which lie behind a crack face (not ahead of the crack front), are recovered "
      "stress-free (stretch, be_bar and plastic deformation gradient set to identity) by the "
      "materials reading this object. Their stored stretch belongs to the smeared crack opening "
      "that the cut crack now carries; recovering it makes the crack faces recoil. Not set: no "
      "relaxation. MOOSE2SIERRA passes the crack cutting threshold by default.");
  params.addParam<std::string>(
      "relax_variable", "d", "Nodal damage variable of the solution file used by relax_threshold");
  params.addParam<std::vector<BoundaryName>>(
      "crack_boundaries",
      "Sidesets of the solution file forming the crack faces for relax_threshold. By default the "
      "external element sides that belong to no sideset (a crack cut before remeshing is not a "
      "sideset of the remeshed file); other free surfaces in no sideset are excluded by their "
      "low damage.");
  return params;
}

SolutionUserObjectQP::SolutionUserObjectQP(const InputParameters & parameters)
  : SolutionUserObject(parameters),
    _tensor_materials(getParam<std::vector<MaterialName>>("tensor_materials")),
    _materials(getParam<std::vector<MaterialName>>("materials")),
    _element(getParam<MooseEnum>("element").getEnum<QpMapping::Element>()),
    _free_after_initial(getParam<bool>("free_after_initial")),
    _payload_freed(false),
    _relax(isParamValid("relax_threshold")),
    _relax_threshold(_relax ? getParam<Real>("relax_threshold") : 1),
    _relax_variable(getParam<std::string>("relax_variable"))
{
  if (!_relax && isParamSetByUser("crack_boundaries"))
    paramError("crack_boundaries", "Only used with 'relax_threshold'.");
  // The damage variable deciding the relaxation is read like the other variables
  if (_relax && std::find(_system_variables.begin(), _system_variables.end(), _relax_variable) ==
                    _system_variables.end())
    _system_variables.push_back(_relax_variable);

  _lookup = QpMapping::getLookup(_element, _qpnum, /*reversed=*/true);

  unsigned int qp_max = _qpnum;
  int dim = 3;

  auto formatQP = [qp_max](unsigned int qp)
  {
    if (qp_max < 10)
      return std::to_string(qp); // Single digit
    else
      return (qp < 10) ? "0" + std::to_string(qp) : std::to_string(qp); // Two digits
  };

  if (!_tensor_materials.empty())
  {
    std::vector<std::string> conv = {"x", "y", "z"};
    for (unsigned int i = 0; i < _tensor_materials.size(); i++)
    {
      std::string matname = _tensor_materials[i];
      if (_tensor_materials[i].size() > 26)
        matname.erase(26, _tensor_materials[i].size() - 1);
      for (unsigned int qp = 1; qp <= qp_max; qp++)
      {
        unsigned int qp_sel = QpMapping::getQP(qp, _lookup);
        for (int j = 0; j < dim; j++)
          for (int k = 0; k < dim; k++)
          {
            _system_variables.push_back(matname + "_" + conv[j] + conv[k] + "_" + formatQP(qp_sel));
          }
      }
    }
  }
  if (!_materials.empty())
  {
    for (unsigned int i = 0; i < _materials.size(); i++)
    {
      // Assuming 8 QPs
      // Starting at 1 because qps start at one in Sierra
      for (unsigned int qp = 1; qp <= qp_max; qp++)
      {
        unsigned int qp_sel = QpMapping::getQP(qp, _lookup);
        _system_variables.push_back(_materials[i] + "_" + formatQP(qp_sel));
      }
    }
  }
}

void
SolutionUserObjectQP::initialSetup()
{
  SolutionUserObject::initialSetup();
  if (_relax && !_relax_locator)
    findRelaxedElements();
}

void
SolutionUserObjectQP::findRelaxedElements()
{
  // The solution file's mesh is replicated, so every processor finds the same elements, also when
  // the restart mesh is distributed.
  const auto & boundary_info = _mesh->get_boundary_info();
  std::set<boundary_id_type> crack_ids;
  const bool by_sideset = isParamValid("crack_boundaries");
  if (by_sideset)
    for (const auto & name : getParam<std::vector<BoundaryName>>("crack_boundaries"))
    {
      const auto id = MooseMeshUtils::getBoundaryID(name, *_mesh);
      if (id == Moose::INVALID_BOUNDARY_ID)
        paramError("crack_boundaries", "Sideset '", name, "' is not in '", _mesh_file, "'.");
      crack_ids.insert(id);
    }

  // Damage of an element: the relax variable at its centroid
  std::unordered_map<dof_id_type, Real> damage;
  const auto element_damage = [&](const Elem & elem)
  {
    const auto it = damage.find(elem.id());
    if (it != damage.end())
      return it->second;
    const Real value = pointValue(_t, elem.vertex_average(), _relax_variable, nullptr);
    damage.emplace(elem.id(), value);
    return value;
  };

  // Crack faces: external sides (in the crack sidesets, or in no sideset) of damaged elements
  struct CrackFace
  {
    Point vertex;
    RealVectorValue edge1, edge2, normal; // normal points out of the material
  };
  std::vector<CrackFace> faces;
  std::vector<const Elem *> seeds;
  std::vector<boundary_id_type> side_ids;
  for (const auto * elem : _mesh->active_element_ptr_range())
    for (const auto s : elem->side_index_range())
    {
      if (elem->neighbor_ptr(s))
        continue;
      boundary_info.boundary_ids(elem, s, side_ids);
      bool crack = side_ids.empty();
      if (by_sideset)
      {
        crack = false;
        for (const auto id : side_ids)
          crack = crack || crack_ids.count(id);
      }
      if (!crack || element_damage(*elem) < _relax_threshold)
        continue;
      const auto side = elem->side_ptr(s);
      const RealVectorValue edge1 = side->point(1) - side->point(0),
                            edge2 = side->point(2) - side->point(0);
      CrackFace face{side->point(0), edge1, edge2, edge1.cross(edge2).unit()};
      if (face.normal * (elem->vertex_average() - face.vertex) > 0)
        face.normal *= -1;
      faces.push_back(face);
      seeds.push_back(elem);
    }

  // A point is behind a crack face when it projects into the face (with some slack for points
  // behind the edges between faces) from the material side.
  const auto behind_crack = [&](const Point & p)
  {
    const Real slack = 0.1;
    for (const auto & face : faces)
    {
      const RealVectorValue r = p - face.vertex;
      if (r * face.normal >= 0)
        continue;
      // Barycentric coordinates of the projection
      const Real a11 = face.edge1 * face.edge1, a12 = face.edge1 * face.edge2,
                 a22 = face.edge2 * face.edge2, b1 = r * face.edge1, b2 = r * face.edge2;
      const Real det = a11 * a22 - a12 * a12;
      const Real u = (a22 * b1 - a12 * b2) / det, v = (a11 * b2 - a12 * b1) / det;
      if (u >= -slack && v >= -slack && 1 - u - v >= -slack)
        return true;
    }
    return false;
  };

  // Elements on the crack faces, and damaged elements connected to them behind the crack faces
  std::vector<const Elem *> front(seeds.begin(), seeds.end());
  for (const auto * elem : seeds)
    _relaxed_elems.insert(elem->id());
  while (!front.empty())
  {
    const auto * elem = front.back();
    front.pop_back();
    for (const auto * neighbor : elem->neighbor_ptr_range())
      if (neighbor && !_relaxed_elems.count(neighbor->id()) &&
          element_damage(*neighbor) >= _relax_threshold &&
          behind_crack(neighbor->vertex_average()))
      {
        _relaxed_elems.insert(neighbor->id());
        front.push_back(neighbor);
      }
  }

  _relax_locator = _mesh->sub_point_locator();
  _relax_locator->enable_out_of_mesh_mode();

  _console << name() << ": relaxing the recovered state of " << _relaxed_elems.size()
           << " elements behind " << faces.size() << " crack faces (" << _relax_variable
           << " >= " << _relax_threshold << ")." << std::endl;
}

bool
SolutionUserObjectQP::relaxRecovery(const Point & p) const
{
  if (!_relax)
    return false;
  if (_payload_freed)
    mooseError(name(), ": relaxRecovery() was called after 'free_after_initial' released the solution.");
  if (_have_last && p.absolute_fuzzy_equals(_last_query))
    return _last_result;
  const Elem * elem = (*_relax_locator)(p);
  _last_query = p;
  _last_result = elem && _relaxed_elems.count(elem->id());
  _have_last = true;
  return _last_result;
}

bool
SolutionUserObjectQP::relaxRecovery(const SolutionUserObject * solution, const Elem & elem)
{
  const auto * qp_solution = dynamic_cast<const SolutionUserObjectQP *>(solution);
  return qp_solution && qp_solution->relaxRecovery(elem.true_centroid());
}

void
SolutionUserObjectQP::timestepSetup()
{
  // Keep the base behavior (cache invalidation, optional time interpolation).
  SolutionUserObject::timestepSetup();

  if (!_free_after_initial || _payload_freed)
    return;

  // Recovery consumers (e.g. ComputeDeformationGradient::initStatefulProperties,
  // LargeDeformationJ2Plasticity*::initQpStatefulProperties) query this object only
  // while stateful material properties are initialized, which completes during
  // FEProblemBase::initialSetup() -- before the first timestepSetup(). By the time we
  // reach here the loaded payload is no longer needed, so release it.
  if (_interpolate_times)
    mooseError("SolutionUserObjectQP '",
               name(),
               "': 'free_after_initial' is incompatible with ExodusII time interpolation, "
               "which needs the loaded solution every step.");

  // Reset in dependency order: the MeshFunctions reference the serialized solution
  // vectors, the systems (owned by the EquationSystems), and the source mesh.
  _relax_locator.reset();
  _mesh_function.reset();
  _mesh_function2.reset();
  _serialized_solution.reset();
  _serialized_solution2.reset();
  _es.reset();
  _es2.reset();
  _system = nullptr;  // raw pointer owned by _es
  _system2 = nullptr; // raw pointer owned by _es2
  _mesh.reset();
  _exodusII_io.reset();
  _nemesis_io.reset();

  _payload_freed = true;

#if defined(__GLIBC__)
  // Freeing the unique_ptrs returns the memory to glibc's allocator, but glibc keeps
  // most of it in its arenas rather than handing it back to the OS. Force a trim so the
  // resident set actually shrinks (the point of free_after_initial).
  malloc_trim(0);
#endif
}
