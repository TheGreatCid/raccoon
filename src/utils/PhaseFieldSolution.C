#include "PhaseFieldSolution.h"
#include "Conversion.h"
#include "MooseObject.h"
#include "MooseEnum.h"

#include "libmesh/elem.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"

InputParameters
PhaseFieldSolution::validParams()
{
  InputParameters params = emptyInputParameters();
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
  return params;
}

PhaseFieldSolution::PhaseFieldSolution(const MooseObject & owner,
                                       const std::vector<std::string> & displacements)
  : _owner(owner), _displacements(displacements), _mesh(owner.comm()), _exodus(_mesh)
{
  const auto & file = owner.getParam<FileName>("solution_file");
  const auto & variable = owner.getParam<std::string>("variable");

  _exodus.read(file);
  _mesh.allow_renumbering(false);
  _mesh.prepare_for_use();

  const auto & nodal_names = _exodus.get_nodal_var_names();
  const auto require = [&](const std::string & param, const std::string & name)
  {
    if (std::find(nodal_names.begin(), nodal_names.end(), name) == nodal_names.end())
      owner.paramError(param,
                       "Nodal variable '",
                       name,
                       "' was not found in '",
                       file,
                       "'. Available: ",
                       Moose::stringify(nodal_names));
  };
  require("variable", variable);
  for (const auto & disp : _displacements)
    require("displacements", disp);

  const int n_steps = _exodus.get_num_time_steps();
  if (n_steps == 0)
    owner.paramError("solution_file", "The file contains no time steps.");
  int timestep = n_steps;
  const auto & timestep_string = owner.getParam<std::string>("timestep");
  if (timestep_string != "LATEST")
  {
    std::istringstream ss(timestep_string);
    if (!((ss >> timestep) && ss.eof()) || timestep < 1 || timestep > n_steps)
      owner.paramError("timestep", "Expected 'LATEST' or an integer from 1 to ", n_steps, ".");
  }

  Order order = FIRST;
  const auto & order_enum = owner.getParam<MooseEnum>("variable_order");
  if (order_enum == "AUTO")
  {
    for (const auto & elem : _mesh.active_element_ptr_range())
      order = std::max(order, elem->default_order());
  }
  else
    order = Utility::string_to_enum<Order>(order_enum);

  _es = std::make_unique<libMesh::EquationSystems>(_mesh);
  auto & system = _es->add_system<libMesh::ExplicitSystem>("phase_field");
  _var_num = system.add_variable(variable, order, LAGRANGE);
  for (const auto & disp : _displacements)
    _disp_nums.push_back(system.add_variable(disp, order, LAGRANGE));
  _es->init();
  _exodus.copy_nodal_solution(system, variable, variable, timestep);
  for (const auto & disp : _displacements)
    _exodus.copy_nodal_solution(system, disp, disp, timestep);

  _serialized = NumericVector<Number>::build(owner.comm());
  _serialized->init(system.n_dofs(), false, SERIAL);
  system.solution->localize(*_serialized);

  _phase_field = std::make_unique<libMesh::MeshFunction>(
      *_es, *_serialized, system.get_dof_map(), _var_num);
  _phase_field->init();
  // Points outside the solution mesh read as undamaged.
  _phase_field->enable_out_of_mesh_mode(Number(0));
}

void
PhaseFieldSolution::moveToDeformed(libMesh::MeshBase & mesh)
{
  mooseAssert(!_disp_nums.empty(), "No displacement variables were read");
  auto & system = _es->get_system("phase_field");

  // Move the input mesh. Its nodes are located in the undeformed solution mesh, so it need not
  // share the file's node numbering.
  {
    libMesh::MeshFunction displacement(*_es, *_serialized, system.get_dof_map(), _disp_nums);
    displacement.init();
    DenseVector<Number> outside(_disp_nums.size(), std::numeric_limits<Real>::max());
    displacement.enable_out_of_mesh_mode(outside);

    std::size_t n_outside = 0;
    DenseVector<Number> u;
    for (auto & node : mesh.node_ptr_range())
    {
      displacement(*node, 0, u);
      if (u(0) == std::numeric_limits<Real>::max())
      {
        ++n_outside;
        continue;
      }
      for (const auto i : index_range(_disp_nums))
        (*node)(i) += u(i);
    }
    if (n_outside)
      _owner.mooseError(n_outside,
                        " input mesh nodes lie outside the undeformed mesh of '",
                        _owner.getParam<FileName>("solution_file"),
                        "', so their displacement is unknown. The input mesh must be in the same "
                        "undeformed configuration as the solution file.");
  }

  // Move the solution mesh by its own nodal displacements, then rebuild the point location.
  const auto sys_num = system.number();
  for (auto & node : _mesh.node_ptr_range())
    for (const auto i : index_range(_disp_nums))
      if (node->n_comp(sys_num, _disp_nums[i]))
        (*node)(i) += (*_serialized)(node->dof_number(sys_num, _disp_nums[i], 0));
  _mesh.clear_point_locator();

  _phase_field = std::make_unique<libMesh::MeshFunction>(
      *_es, *_serialized, system.get_dof_map(), _var_num);
  _phase_field->init();
  _phase_field->enable_out_of_mesh_mode(Number(0));
}

std::pair<Real, Real>
PhaseFieldSolution::findRidge(const Point & p, const RealVectorValue & unit_normal, const Real reach)
{
  // Sample densely enough to resolve the quadratic variation of d within each element crossed.
  const int n = 20;
  Real peak = -std::numeric_limits<Real>::max();
  Real offset = 0;
  for (int k = -n; k <= n; ++k)
  {
    const Real s = reach * k / n;
    const Real value = (*_phase_field)(p + s * unit_normal);
    // Prefer the sample closest to p among equal values, e.g. on a d = 1 plateau
    if (value > peak || (value == peak && std::abs(s) < std::abs(offset)))
    {
      peak = value;
      offset = s;
    }
  }
  return {offset, peak};
}

Real
PhaseFieldSolution::ridgeReach(const std::vector<const Elem *> & elems, const Real distance)
{
  if (distance > 0)
    return distance;
  Real reach = 0;
  for (const auto * elem : elems)
    reach = std::max(reach, 2 * elem->hmax());
  return reach;
}
