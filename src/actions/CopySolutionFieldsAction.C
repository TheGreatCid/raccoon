#include "CopySolutionFieldsAction.h"
#include "Conversion.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseMesh.h"

#include "libmesh/exodusII_io_helper.h"

registerMooseAction("raccoonApp", CopySolutionFieldsAction, "add_aux_variable");
registerMooseAction("raccoonApp", CopySolutionFieldsAction, "add_user_object");
registerMooseAction("raccoonApp", CopySolutionFieldsAction, "add_ic");

InputParameters
CopySolutionFieldsAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription(
      "Copies all nodal and elemental variables of an Exodus file onto the mesh as AuxVariables "
      "set by SolutionIC, reading the variable names from the file.");
  params.addRequiredParam<FileName>("file", "The Exodus file to copy the variables from");
  params.addParam<std::string>(
      "timestep", "LATEST", "The Exodus time step (1-based) to copy, or 'LATEST'");
  params.addParam<std::vector<std::string>>(
      "variables", "Copy only these variables. By default every nodal and elemental variable.");
  params.addParam<std::vector<std::string>>("exclude", {}, "Variables not to copy");
  params.addParam<MeshGeneratorName>(
      "moved_nodes_from",
      "A mesh generator that moved nodes and recorded it (e.g. PlanarCrackGenerator). Values are "
      "then looked up at the nodes' positions before the move, so moved nodes keep their original "
      "values.");
  return params;
}

CopySolutionFieldsAction::CopySolutionFieldsAction(const InputParameters & params)
  : Action(params), _file(getParam<FileName>("file")), _header_read(false)
{
}

std::string
CopySolutionFieldsAction::auxVariableName(const std::string & file_variable) const
{
  return file_variable;
}

void
CopySolutionFieldsAction::readHeader()
{
  if (_header_read)
    return;
  _header_read = true;

  // Read on every processor so each one gets the names without communication.
  libMesh::ExodusII_IO_Helper exodus(*this, false, /*run_only_on_proc0=*/false);
  exodus.open(_file.c_str(), /*read_only=*/true);
  exodus.read_and_store_header_info();
  exodus.read_block_info();
  for (int i = 0; i < exodus.num_elem_blk; ++i)
  {
    const auto block_name = exodus.get_block_name(i);
    _file_blocks.push_back(block_name.empty() ? std::to_string(exodus.get_block_id(i))
                                              : block_name);
  }
  exodus.read_var_names(libMesh::ExodusII_IO_Helper::NODAL);
  exodus.read_var_names(libMesh::ExodusII_IO_Helper::ELEMENTAL);
  const auto nodal = exodus.nodal_var_names;
  const auto elemental = exodus.elem_var_names;
  exodus.close();

  const auto contains = [](const std::vector<std::string> & names, const std::string & name)
  { return std::find(names.begin(), names.end(), name) != names.end(); };

  std::set<std::string> selected;
  if (isParamValid("variables"))
    for (const auto & var : getParam<std::vector<std::string>>("variables"))
    {
      if (!contains(nodal, var) && !contains(elemental, var))
        paramError("variables",
                   "Variable '",
                   var,
                   "' is not a nodal or elemental variable of '",
                   _file,
                   "'.");
      selected.insert(var);
    }

  const auto & exclude = getParam<std::vector<std::string>>("exclude");
  for (const auto & var : exclude)
    if (!contains(nodal, var) && !contains(elemental, var))
      paramError(
          "exclude", "Variable '", var, "' is not a nodal or elemental variable of '", _file, "'.");

  const auto keep = [&](const std::string & var)
  { return (selected.empty() || selected.count(var)) && !contains(exclude, var); };

  for (const auto & var : nodal)
    if (keep(var))
      _nodal_variables.push_back(var);
  for (const auto & var : elemental)
    if (keep(var))
    {
      // A name used for both a nodal and an elemental variable cannot become two AuxVariables.
      if (contains(_nodal_variables, var))
        mooseError(name(),
                   ": '",
                   var,
                   "' is both a nodal and an elemental variable in '",
                   _file,
                   "'; add it to 'exclude'.");
      _elemental_variables.push_back(var);
    }

  if (_nodal_variables.empty() && _elemental_variables.empty())
    mooseWarning(name(), ": no variables to copy from '", _file, "'.");
}

void
CopySolutionFieldsAction::act()
{
  readHeader();
  const std::string solution_name = name() + "_solution";

  if (_current_task == "add_aux_variable")
  {
    const std::string nodal_order = _mesh->hasSecondOrderElements() ? "SECOND" : "FIRST";
    for (const auto & variable : _nodal_variables)
    {
      auto params = _factory.getValidParams("MooseVariable");
      params.set<MooseEnum>("family") = "LAGRANGE";
      params.set<MooseEnum>("order") = nodal_order;
      _problem->addAuxVariable("MooseVariable", auxVariableName(variable), params);
    }
    for (const auto & variable : _elemental_variables)
    {
      auto params = _factory.getValidParams("MooseVariable");
      params.set<MooseEnum>("family") = "MONOMIAL";
      params.set<MooseEnum>("order") = "CONSTANT";
      _problem->addAuxVariable("MooseVariable", auxVariableName(variable), params);
    }
  }

  else if (_current_task == "add_user_object")
  {
    if (_nodal_variables.empty() && _elemental_variables.empty())
      return;
    auto params = _factory.getValidParams("SolutionUserObject");
    params.set<MeshFileName>("mesh") = _file;
    params.set<std::string>("timestep") = getParam<std::string>("timestep");
    // Second-order meshes need the mid-edge values, not a linear interpolation of the vertices.
    params.set<MooseEnum>("nodal_variable_order") =
        _mesh->hasSecondOrderElements() ? "SECOND" : "FIRST";
    auto & variables = params.set<std::vector<std::string>>("system_variables");
    variables.assign(_nodal_variables.begin(), _nodal_variables.end());
    variables.insert(variables.end(), _elemental_variables.begin(), _elemental_variables.end());
    _problem->addUserObject("SolutionUserObject", solution_name, params);
  }

  else if (_current_task == "add_ic")
  {
    const bool moved = isParamValid("moved_nodes_from");
    const std::string ic_type = moved ? "OriginalPositionSolutionIC" : "SolutionIC";
    std::vector<Point> moved_positions, moved_offsets;
    if (moved)
    {
      const auto & generator = getParam<MeshGeneratorName>("moved_nodes_from");
      for (const auto & property : {"moved_node_positions", "moved_node_offsets"})
        if (!hasMeshProperty<std::vector<Point>>(property, generator))
          paramError("moved_nodes_from",
                     "Mesh generator '",
                     generator,
                     "' did not record moved nodes ('",
                     property,
                     "').");
      moved_positions = getMeshProperty<std::vector<Point>>("moved_node_positions", generator);
      moved_offsets = getMeshProperty<std::vector<Point>>("moved_node_offsets", generator);
    }

    const auto add_ic = [&](const std::string & variable)
    {
      auto params = _factory.getValidParams(ic_type);
      if (moved)
      {
        params.set<std::vector<Point>>("moved_node_positions") = moved_positions;
        params.set<std::vector<Point>>("moved_node_offsets") = moved_offsets;
      }
      params.set<VariableName>("variable") = auxVariableName(variable);
      params.set<UserObjectName>("solution_uo") = solution_name;
      params.set<VariableName>("from_variable") = variable;
      // The mesh may have subdomains the file does not (e.g. crack blocks), so read from all of
      // the file's blocks instead of matching subdomain names.
      params.set<std::vector<SubdomainName>>("from_subdomains") = _file_blocks;
      _problem->addInitialCondition(ic_type, name() + "_" + variable, params);
    };
    for (const auto & variable : _nodal_variables)
      add_ic(variable);
    for (const auto & variable : _elemental_variables)
      add_ic(variable);
  }
}
