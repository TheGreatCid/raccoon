//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

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
  params.addParam<std::vector<MaterialName>>(
      "tensor_materials", {}, "Tensor properties to read at the quadrature points");
  params.addParam<std::vector<MaterialName>>(
      "materials", {}, "Scalar properties to read at the quadrature points");
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
  return params;
}

SolutionUserObjectQP::SolutionUserObjectQP(const InputParameters & parameters)
  : SolutionUserObject(parameters),
    _tensor_materials(getParam<std::vector<MaterialName>>("tensor_materials")),
    _materials(getParam<std::vector<MaterialName>>("materials")),
    _element(getParam<MooseEnum>("element").getEnum<QpMapping::Element>()),
    _free_after_initial(getParam<bool>("free_after_initial")),
    _payload_freed(false)
{
  _lookup = QpMapping::getLookup(_element, _qpnum, /*reversed=*/true);

  // Read the variables of the listed properties at every quadrature point
  for (const auto qp : make_range(_qpnum))
  {
    for (const auto & material : _materials)
      _system_variables.push_back(qpVariableName(material, qp));
    for (const auto & material : _tensor_materials)
      for (const auto i : make_range(3u))
        for (const auto j : make_range(3u))
          _system_variables.push_back(qpVariableName(material, i, j, qp));
  }
}

std::string
SolutionUserObjectQP::qpVariableName(const std::string & property, const unsigned int qp) const
{
  return QpMapping::variableName(property, QpMapping::getQP(qp + 1, _lookup), _qpnum);
}

std::string
SolutionUserObjectQP::qpVariableName(const std::string & property,
                                     const unsigned int i,
                                     const unsigned int j,
                                     const unsigned int qp) const
{
  return QpMapping::variableName(property, i, j, QpMapping::getQP(qp + 1, _lookup), _qpnum);
}

Real
SolutionUserObjectQP::qpValue(const Real t,
                              const Elem & elem,
                              const std::string & property,
                              const unsigned int qp) const
{
  return pointValue(t, elem.true_centroid(), qpVariableName(property, qp));
}

RankTwoTensor
SolutionUserObjectQP::qpTensor(const Real t,
                               const Elem & elem,
                               const std::string & property,
                               const unsigned int qp) const
{
  const Point centroid = elem.true_centroid();
  RankTwoTensor tensor;
  for (const auto i : make_range(3u))
    for (const auto j : make_range(3u))
      tensor(i, j) = pointValue(t, centroid, qpVariableName(property, i, j, qp));
  return tensor;
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
