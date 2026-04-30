//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "MidEdgeNodalAverager.h"
#include "AuxiliarySystem.h"
#include "MooseMesh.h"
#include "MooseVariableFE.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"

registerMooseObject("raccoonApp", MidEdgeNodalAverager);

InputParameters
MidEdgeNodalAverager::validParams()
{
  InputParameters p = GeneralUserObject::validParams();
  p.addRequiredParam<std::vector<AuxVariableName>>(
      "variables",
      "Second-order nodal AuxVariables whose mid-edge values should be replaced "
      "with the average of their two end-node values.");
  ExecFlagEnum & exec = p.set<ExecFlagEnum>("execute_on", true);
  exec = EXEC_INITIAL;
  return p;
}

MidEdgeNodalAverager::MidEdgeNodalAverager(const InputParameters & params)
  : GeneralUserObject(params),
    _var_names(getParam<std::vector<AuxVariableName>>("variables"))
{
}

void
MidEdgeNodalAverager::execute()
{
  auto & aux_sys = _fe_problem.getAuxiliarySystem();
  NumericVector<Number> & sol = aux_sys.solution();
  const NumericVector<Number> & sol_ghost = *aux_sys.currentSolution();
  const unsigned int sys_num = aux_sys.number();
  MooseMesh & mesh = _subproblem.mesh();

  for (const auto & vname : _var_names)
  {
    auto & var = aux_sys.getFieldVariable<Real>(0, vname);
    const unsigned int var_num = var.number();

    for (const auto & elem : *mesh.getActiveLocalElementRange())
    {
      const unsigned int n_corners = elem->n_vertices();
      for (unsigned int e = 0; e < elem->n_edges(); ++e)
      {
        const auto edge_nodes = elem->nodes_on_edge(e);

        unsigned int n_c = 0;
        unsigned int mid_local = libMesh::invalid_uint;
        unsigned int c_local[2] = {0, 0};
        for (auto l : edge_nodes)
        {
          if (l < n_corners)
          {
            if (n_c < 2)
              c_local[n_c++] = l;
          }
          else
            mid_local = l;
        }
        if (n_c != 2 || mid_local == libMesh::invalid_uint)
          continue;

        const Node & nm = elem->node_ref(mid_local);
        if (nm.processor_id() != processor_id())
          continue;

        const Node & n0 = elem->node_ref(c_local[0]);
        const Node & n1 = elem->node_ref(c_local[1]);

        const dof_id_type d0 = n0.dof_number(sys_num, var_num, 0);
        const dof_id_type d1 = n1.dof_number(sys_num, var_num, 0);
        const dof_id_type dm = nm.dof_number(sys_num, var_num, 0);

        sol.set(dm, 0.5 * (sol_ghost(d0) + sol_ghost(d1)));
      }
    }
  }

  sol.close();
  aux_sys.update();
}
