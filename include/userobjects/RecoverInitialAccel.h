//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ElementUserObject.h"

#include <unordered_map>
#include <set>

/**
 * Seeds a consistent initial (Newmark/HHT) acceleration for a dynamic restart.
 *
 * A recovered state generally is not in discrete equilibrium at the restart time: reading the
 * acceleration from the dump leaves M*a0 + F_int(u0) - F_ext(t0) != 0, and the first time step
 * corrects that imbalance -- which the second time-derivative amplifies into an acceleration spike
 * (and rings transverse modes). This user object instead computes the acceleration that the
 * recovered forces actually produce, using a lumped mass so no linear solve is needed:
 *
 *     a_i = - F_int,i / m_i
 *     F_int,i = \int (grad(phi_i) . sigma) dV     (internal force from `stress`)
 *     m_i     = \int rho phi_i dV                 (lumped mass from `density`)
 *
 * Integration is on the displaced mesh (Cauchy stress). It runs once on INITIAL and overwrites the
 * acceleration aux variables on the free (non-`skip_boundaries`) nodes so the first HHT/Newmark
 * step starts in balance.
 */
class RecoverInitialAccel : public ElementUserObject
{
public:
  static InputParameters validParams();

  RecoverInitialAccel(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & uo) override;
  virtual void finalize() override;

protected:
  /// Displacement variables (FE type / geometry / node layout)
  const unsigned int _ndisp;
  std::vector<MooseVariable *> _disp_var;
  /// Acceleration aux variables to seed (one per displacement component)
  std::vector<MooseVariable *> _accel_var;

  /// Cauchy stress used to build the internal force
  const ADMaterialProperty<RankTwoTensor> & _stress;
  /// Mass density used to build the lumped mass
  const ADMaterialProperty<Real> & _density;

  unsigned int _dim;

  /// Per-component sets of nodes whose acceleration is left untouched (the component is constrained
  /// / prescribed there, so its acceleration is imposed by the boundary condition).
  std::vector<std::set<dof_id_type>> _skip_nodes;

  /// Collect the nodes lying on the given boundaries (handles both nodesets and sidesets).
  std::set<dof_id_type> nodesOnBoundaries(const std::vector<BoundaryName> & names) const;

  /// Thread/rank-local accumulators keyed by acceleration dof
  std::unordered_map<dof_id_type, Real> _fint;
  std::unordered_map<dof_id_type, Real> _mass;
};
