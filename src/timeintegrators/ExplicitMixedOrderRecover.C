//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ExplicitMixedOrderRecover.h"
#include "SolutionUserObject.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/system.h"
#include "libmesh/node.h"

registerMooseObject("raccoonApp", ExplicitMixedOrderRecover);

InputParameters
ExplicitMixedOrderRecover::validParams()
{
  InputParameters params = ExplicitMixedOrder::validParams();
  params.addClassDescription(
      "Recovery-aware ExplicitMixedOrder.  After the base class init() runs "
      "(which zeros solutionUDot via the forward Euler estimate), reads "
      "recovered nodal velocity / acceleration values from a SolutionUserObject "
      "pointing at the reference dump and writes them directly into the "
      "integrator's solutionUDot / solutionUDotDot vectors at INITIAL.  Lets "
      "the central-difference solve pick up at the reference's state instead "
      "of starting from zero velocity.");
  params.addRequiredParam<UserObjectName>(
      "solution",
      "The SolutionUserObject reading the reference's recovery exodus dump.");
  params.addRequiredParam<std::vector<VariableName>>(
      "disp_vars",
      "Nonlinear displacement variables on the restart side (target DOFs for "
      "the recovered vel / accel writes).");
  params.addRequiredParam<std::vector<std::string>>(
      "vel_vars",
      "Names of nodal velocity AuxVariables in the dump, one per disp_vars entry.");
  params.addRequiredParam<std::vector<std::string>>(
      "accel_vars",
      "Names of nodal acceleration AuxVariables in the dump, one per disp_vars entry.");
  return params;
}

ExplicitMixedOrderRecover::ExplicitMixedOrderRecover(const InputParameters & parameters)
  : ExplicitMixedOrder(parameters),
    _sol_uo(nullptr),
    _disp_var_names(getParam<std::vector<VariableName>>("disp_vars")),
    _vel_var_names(getParam<std::vector<std::string>>("vel_vars")),
    _accel_var_names(getParam<std::vector<std::string>>("accel_vars")),
    _first_solve_pending(true)
{
  if (_disp_var_names.size() != _vel_var_names.size() ||
      _disp_var_names.size() != _accel_var_names.size())
    mooseError("disp_vars (",
               _disp_var_names.size(),
               "), vel_vars (",
               _vel_var_names.size(),
               "), accel_vars (",
               _accel_var_names.size(),
               ") must all have the same length.");
}

void
ExplicitMixedOrderRecover::init()
{
  // Base class runs first: sets up first/second-order DOF index lists and
  // calls computeICs(), which zeros solutionUDot via the (u - u_old)/dt
  // estimate.  We override that zero immediately afterwards.
  ExplicitMixedOrder::init();

  _sol_uo = &_fe_problem.getUserObject<SolutionUserObject>(
      getParam<UserObjectName>("solution"));

  recoverVelAccel();
}

void
ExplicitMixedOrderRecover::preSolve()
{
  ExplicitMixedOrder::preSolve();

  // On the very first restart step, _dt_old is whatever the executioner
  // initialized it to (typically 0).  ExplicitMixedOrder's central-difference
  // velocity update scales the new acceleration by (_dt + _dt_old)/2:
  //
  //   vel_{n+1} = vel_n + a_n * (dt + dt_old)/2
  //
  // If _dt_old = 0, only half the accel contribution is applied -- the
  // restart's first step underintegrates the velocity vs. the equivalent
  // step in a continuous run (where dt_old = previous dt).  Set them equal
  // here so the first restart step matches what the reference's continued
  // run would compute.  Subsequent steps use the executioner's normal
  // dt_old advancement.
  if (_first_solve_pending)
  {
    _dt_old = _dt;
    _first_solve_pending = false;
    _console << "[ExplicitMixedOrderRecover] first restart step: forced "
                "_dt_old = _dt = "
             << _dt << " so the CD velocity update uses the full "
                       "acceleration contribution."
             << std::endl;
  }
}

void
ExplicitMixedOrderRecover::recoverVelAccel()
{
  // Current-step time derivative vectors.
  auto * vel = _sys.solutionUDot();
  auto * accel = _sys.solutionUDotDot();
  if (!vel || !accel)
    mooseError("ExplicitMixedOrderRecover: solutionUDot / solutionUDotDot is "
               "null.  Base ExplicitMixedOrder constructor should request "
               "them via setUDotRequested(true) / setUDotDotRequested(true).");

  // Old-step time derivative vectors.  ExplicitDirichletBCBase reads
  // nodalValueDotOld() in its first-step residual, which pulls from
  // solutionUDotOld().  If we don't seed them, the BC's _u_dot_old = 0 at
  // the first restart step and it computes a residual as if the body were
  // at rest -- causing a large transient kick on the first step.
  auto * vel_old = _sys.solutionUDotOld();
  auto * accel_old = _sys.solutionUDotDotOld();

  auto & lm_sys = _sys.system();
  const auto & mesh = lm_sys.get_mesh();
  const auto sys_num = lm_sys.number();

  unsigned int n_nodes_written = 0;

  for (std::size_t i = 0; i < _disp_var_names.size(); ++i)
  {
    const auto & disp_name = _disp_var_names[i];
    if (!lm_sys.has_variable(disp_name))
      mooseError("ExplicitMixedOrderRecover: nonlinear system has no "
                 "variable named '",
                 disp_name,
                 "'.");
    const auto var_num = lm_sys.variable_number(disp_name);

    const auto & vel_name = _vel_var_names[i];
    const auto & accel_name = _accel_var_names[i];

    // Iterate local nodes (owned by this rank).  For each node that carries
    // a DOF for the disp variable, write the recovered vel / accel into the
    // integrator vectors at that DOF.
    for (const auto * node : as_range(mesh.local_nodes_begin(), mesh.local_nodes_end()))
    {
      const auto n_comp = node->n_comp(sys_num, var_num);
      if (n_comp == 0)
        continue;

      // disp variables are scalar -> a single component per node (comp 0).
      const auto dof = node->dof_number(sys_num, var_num, 0);
      if (dof == libMesh::DofObject::invalid_id)
        continue;

      const Real v = _sol_uo->directValue(node, vel_name);
      const Real a = _sol_uo->directValue(node, accel_name);
      vel->set(dof, v);
      accel->set(dof, a);
      // Seed the "old" vectors too so the BC's first-step formula sees the
      // recovered vel/accel as _u_dot_old / _u_dotdot_old.
      if (vel_old)
        vel_old->set(dof, v);
      if (accel_old)
        accel_old->set(dof, a);
      ++n_nodes_written;
    }
  }

  vel->close();
  accel->close();
  if (vel_old)
    vel_old->close();
  if (accel_old)
    accel_old->close();

  _console << "[ExplicitMixedOrderRecover] seeded solutionUDot / "
              "solutionUDotDot (and their Old vectors) from dump at "
           << n_nodes_written << " (local node, disp-component) pairs."
           << std::endl;
}
