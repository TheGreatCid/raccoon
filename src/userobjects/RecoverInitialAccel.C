//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "RecoverInitialAccel.h"

#include "AuxiliarySystem.h"
#include "MooseVariable.h"
#include "MooseMesh.h"

#include "metaphysicl/raw_type.h"
#include "libmesh/mesh_base.h"

registerMooseObject("raccoonApp", RecoverInitialAccel);

InputParameters
RecoverInitialAccel::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params.addClassDescription(
      "Seeds a consistent initial (Newmark/HHT) acceleration for a dynamic restart by solving the "
      "lumped equation of motion a_i = -F_int,i / m_i at the recovered state, where F_int is the "
      "internal force built from `stress` and m is the lumped mass built from `density`. Runs once "
      "on INITIAL so the first time step starts in balance and does not kick the acceleration.");
  params.addRequiredCoupledVar("displacements", "The displacement variables");
  params.addRequiredCoupledVar("accelerations",
                               "The acceleration aux variables to seed (one per displacement).");
  params.addParam<MaterialPropertyName>(
      "stress",
      "stress",
      "Cauchy stress used to build the internal force. Use the constitutive `stress` (what the "
      "residual evaluates) so the seeded acceleration balances the recovered state and the first "
      "step does not kick.");
  params.addParam<MaterialPropertyName>(
      "density",
      "density",
      "Mass density used to build the lumped mass. For a recover run set this to the same density "
      "the inertia kernels use (e.g. adj_density_rec).");
  params.addParam<std::vector<BoundaryName>>(
      "skip_x_boundaries",
      {},
      "Boundaries (nodesets or sidesets) on which disp_x is constrained/prescribed; accel_x there is "
      "left untouched (imposed by its BC).");
  params.addParam<std::vector<BoundaryName>>(
      "skip_y_boundaries", {}, "Boundaries on which disp_y is constrained; accel_y left untouched.");
  params.addParam<std::vector<BoundaryName>>(
      "skip_z_boundaries", {}, "Boundaries on which disp_z is constrained; accel_z left untouched.");

  // Cauchy stress -> integrate on the displaced mesh; seed once at the restart.
  params.set<bool>("use_displaced_mesh") = true;
  ExecFlagEnum & exec = params.set<ExecFlagEnum>("execute_on", true);
  exec = {EXEC_INITIAL};
  return params;
}

RecoverInitialAccel::RecoverInitialAccel(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _ndisp(coupledComponents("displacements")),
    _stress(getADMaterialProperty<RankTwoTensor>("stress")),
    _density(getADMaterialProperty<Real>("density"))
{
  if (coupledComponents("accelerations") != _ndisp)
    paramError("accelerations",
               "The number of acceleration variables must match the number of displacements.");
  for (unsigned int c = 0; c < _ndisp; ++c)
  {
    _disp_var.push_back(getVar("displacements", c));
    _accel_var.push_back(getVar("accelerations", c));
  }
  _dim = _mesh.dimension();
}

std::set<dof_id_type>
RecoverInitialAccel::nodesOnBoundaries(const std::vector<BoundaryName> & names) const
{
  std::set<dof_id_type> out;
  if (names.empty())
    return out;

  std::set<boundary_id_type> ids;
  for (const auto & bn : names)
    ids.insert(_mesh.getBoundaryID(bn));

  const auto & binfo = _mesh.getMesh().get_boundary_info();

  // Nodeset nodes.
  for (const auto & [nid, bid] : binfo.build_node_list())
    if (ids.count(bid))
      out.insert(nid);

  // Sideset nodes (for boundaries defined as sidesets, e.g. a prescribed-displacement face).
  for (const auto & tup : binfo.build_side_list())
  {
    const auto elem_id = std::get<0>(tup);
    const auto side = std::get<1>(tup);
    const auto bid = std::get<2>(tup);
    if (!ids.count(bid))
      continue;
    const Elem * elem = _mesh.getMesh().elem_ptr(elem_id);
    for (const auto n : elem->nodes_on_side(side))
      out.insert(elem->node_id(n));
  }
  return out;
}

void
RecoverInitialAccel::initialSetup()
{
  // Resolve the per-component constrained nodes once.
  _skip_nodes.assign(_ndisp, std::set<dof_id_type>());
  const char * comp_param[3] = {"skip_x_boundaries", "skip_y_boundaries", "skip_z_boundaries"};
  for (unsigned int c = 0; c < _ndisp && c < 3; ++c)
    _skip_nodes[c] = nodesOnBoundaries(getParam<std::vector<BoundaryName>>(comp_param[c]));
}

void
RecoverInitialAccel::initialize()
{
  _fint.clear();
  _mass.clear();
}

void
RecoverInitialAccel::execute()
{
  // Use the assembly's own shape functions/gradients for the displacement variable: they are
  // evaluated at exactly the material's quadrature points on the current (displaced) element, so
  // sigma[qp] and grad(phi_i)[qp] are paired consistently (a standalone FE risks a QP mismatch).
  const VariablePhiValue & phi = _disp_var[0]->phi();
  const VariablePhiGradient & grad_phi = _disp_var[0]->gradPhi();

  const unsigned int nqp = _qrule->n_points();
  const unsigned int nsh = phi.size();
  const unsigned int aux_sys_num = _accel_var[0]->sys().number();

  for (unsigned int i = 0; i < nsh; ++i)
  {
    const Node & node = _current_elem->node_ref(i);
    for (unsigned int c = 0; c < _ndisp; ++c)
    {
      const dof_id_type dof = node.dof_number(aux_sys_num, _accel_var[c]->number(), 0);
      Real fint = 0.0;
      Real mass = 0.0;
      for (unsigned int qp = 0; qp < nqp; ++qp)
      {
        const RankTwoTensor sig = MetaPhysicL::raw_value(_stress[qp]);
        const Real rho = MetaPhysicL::raw_value(_density[qp]);
        const Real w = _JxW[qp] * _coord[qp];
        // Internal nodal force F_int_{i,c} = \int sigma.row(c) . grad(phi_i) dV
        fint += (sig.row(c) * grad_phi[i][qp]) * w;
        // Lumped nodal mass m_{i} = \int rho phi_i dV
        mass += rho * phi[i][qp] * w;
      }
      _fint[dof] += fint;
      _mass[dof] += mass;
    }
  }
}

void
RecoverInitialAccel::threadJoin(const UserObject & uo)
{
  const auto & other = static_cast<const RecoverInitialAccel &>(uo);
  for (const auto & [dof, v] : other._fint)
    _fint[dof] += v;
  for (const auto & [dof, v] : other._mass)
    _mass[dof] += v;
}

void
RecoverInitialAccel::finalize()
{
  auto & aux_sys = _fe_problem.getAuxiliarySystem();
  NumericVector<Number> & sol = aux_sys.solution();

  // Sum contributions across ranks: push the (thread-joined) rank-local contributions into
  // parallel vectors and close with ADD so shared/partition-boundary nodes are summed.
  auto fint_vec = sol.zero_clone();
  auto mass_vec = sol.zero_clone();
  for (const auto & [dof, v] : _fint)
    fint_vec->add(dof, v);
  for (const auto & [dof, v] : _mass)
    mass_vec->add(dof, v);
  fint_vec->close();
  mass_vec->close();

  // Write a_i = -F_int,i / m_i on locally-owned, non-skipped acceleration dofs.
  const unsigned int aux_sys_num = aux_sys.number();
  unsigned long n_written = 0;
  Real max_abs = 0.0;
  for (const auto & node : _mesh.getMesh().local_node_ptr_range())
  {
    for (unsigned int c = 0; c < _ndisp; ++c)
    {
      if (_skip_nodes[c].count(node->id()))
        continue;
      if (node->n_dofs(aux_sys_num, _accel_var[c]->number()) == 0)
        continue;
      const dof_id_type dof = node->dof_number(aux_sys_num, _accel_var[c]->number(), 0);
      const Real m = (*mass_vec)(dof);
      if (m > 0.0)
      {
        const Real a = -(*fint_vec)(dof) / m;
        sol.set(dof, a);
        ++n_written;
        max_abs = std::max(max_abs, std::abs(a));
      }
    }
  }
  sol.close();
  _communicator.sum(n_written);
  _communicator.max(max_abs);
  std::size_t n_skip = 0;
  for (const auto & s : _skip_nodes)
    n_skip += s.size();
  _console << "RecoverInitialAccel '" << name() << "': seeded consistent a0 on " << n_written
           << " free accel dofs (max|a0| = " << max_abs << ", " << n_skip
           << " constrained dofs left to their BCs)." << std::endl;

  // Push the updated acceleration into the ghosted/current solution the kernels read.
  aux_sys.update();
}
