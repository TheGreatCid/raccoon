#pragma once

#include "SolutionUserObject.h"
#include "Qp_Mapping.h"

#include "libmesh/point_locator_base.h"

#include <unordered_set>

class SolutionUserObjectQP : public SolutionUserObject
{
public:
  static InputParameters validParams();
  SolutionUserObjectQP(const InputParameters & parameters);

  virtual void initialSetup() override;

  /// Whether this object relaxes the recovered state of some elements (relax_threshold is set)
  bool relaxing() const { return _relax; }

  /**
   * Whether the recovered state at a point (a restart element's centroid) is relaxed to a
   * stress-free state, because the point lies in a damaged element behind a cut crack face. Only
   * valid while the recovery is read, before a free_after_initial release.
   */
  bool relaxRecovery(const Point & p) const;

  /// relaxRecovery() for any SolutionUserObject: false unless it is a relaxing SolutionUserObjectQP
  static bool relaxRecovery(const SolutionUserObject * solution, const Elem & elem);

protected:
  /// Releases the large resident payload after the one-time initial recovery reads
  virtual void timestepSetup() override;

  /// Finds the elements whose recovered state is relaxed (see relax_threshold)
  void findRelaxedElements();

  std::vector<MaterialName> _tensor_materials;
  std::vector<MaterialName> _materials;

  QpMapping::Element _element;

  unsigned int _qpnum;

  /// If true, release the loaded solution/mesh payload after the one-time initial
  /// recovery reads. Only safe for the recovery pattern (this object is queried
  /// solely during stateful-property initialization); incompatible with time
  /// interpolation. Do not enable when the object is queried after INITIAL.
  const bool _free_after_initial;

  /// Guard so the payload is released only once
  bool _payload_freed;

  /// Whether the recovered state of damaged elements behind the crack faces is relaxed
  const bool _relax;
  /// Elements behind the crack faces are relaxed where the relax variable reaches this
  const Real _relax_threshold;
  /// Nodal damage variable of the solution file deciding which elements are relaxed
  const std::string _relax_variable;

  /// Relaxed elements of the solution file's mesh (ids of _mesh)
  std::unordered_set<dof_id_type> _relaxed_elems;
  /// Locates restart points in the solution file's mesh for relaxRecovery()
  std::unique_ptr<libMesh::PointLocatorBase> _relax_locator;
  /// Last relaxRecovery() query, which consumers repeat for every quadrature point of an element
  mutable Point _last_query;
  mutable bool _last_result = false;
  mutable bool _have_last = false;

private:
  const std::unordered_map<int, int> * _lookup;
};
