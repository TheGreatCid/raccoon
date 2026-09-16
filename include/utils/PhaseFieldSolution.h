#pragma once

#include "InputParameters.h"
#include "MooseTypes.h"

#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/replicated_mesh.h"

class MooseObject;

/**
 * A nodal phase-field variable (and optionally displacements) read from one time step of an Exodus
 * file, evaluable anywhere through a MeshFunction. Used by the mesh generators that turn a
 * phase-field crack into a sharp crack, which run before any UserObject exists.
 */
class PhaseFieldSolution
{
public:
  /// Parameters for the file, variable, time step and interpolation order
  static InputParameters validParams();

  /**
   * Reads the file named by the owner's parameters. Parameter errors are reported through the
   * owner. When displacement variables are given they are read as well (see moveToDeformed).
   */
  PhaseFieldSolution(const MooseObject & owner,
                     const std::vector<std::string> & displacements = {});

  /**
   * Moves the nodes of a mesh in the file's undeformed configuration by the file's displacements,
   * and moves the solution mesh too, so later evaluations are in the deformed configuration.
   * Requires displacements to have been given.
   */
  void moveToDeformed(libMesh::MeshBase & mesh);

  /// The phase field at p (0 outside the solution mesh)
  Real value(const Point & p) { return (*_phase_field)(p); }
  /// The phase-field gradient at p
  RealVectorValue gradient(const Point & p) { return _phase_field->gradient(p); }

  /**
   * The ridge along a line through p: its signed offset from p along unit_normal and the largest
   * phase-field value, searched within reach on each side of p
   */
  std::pair<Real, Real>
  findRidge(const Point & p, const RealVectorValue & unit_normal, Real reach);

  /// Twice the largest size of the elements, or the given distance when it is positive
  static Real ridgeReach(const std::vector<const Elem *> & elems, Real distance);

private:
  const MooseObject & _owner;
  const std::vector<std::string> _displacements;
  libMesh::ReplicatedMesh _mesh;
  libMesh::ExodusII_IO _exodus;
  std::unique_ptr<libMesh::EquationSystems> _es;
  std::unique_ptr<NumericVector<Number>> _serialized;
  std::unique_ptr<libMesh::MeshFunction> _phase_field;
  unsigned int _var_num;
  std::vector<unsigned int> _disp_nums;
};
