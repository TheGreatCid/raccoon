#include "Factory.h"
#include "FEProblem.h"
#include "AddVariableAction.h"
#include "InputParameters.h"
#include "NonlinearSystemBase.h"
#include "Parser.h"
#include "RecoverVariablesAction.h"
#include <vector>

registerMooseAction("raccoonApp", RecoverVariablesAction, "add_aux_variable");
registerMooseAction("raccoonApp", RecoverVariablesAction, "add_aux_kernel");

InputParameters
RecoverVariablesAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addParam<std::vector<MaterialName>>("tensor_materials", "materials to output qps on");
  params.addParam<std::vector<MaterialName>>("materials", "materials to output qps on");
  params.addParam<MaterialName>("def_grad_name", "name of deformation gradient");
  params.addParam<unsigned int>("recover_num", 0, "number of recovery");
  return params;
}

RecoverVariablesAction::RecoverVariablesAction(const InputParameters & params)
  : Action(params),
    _tensor_materials(getParam<std::vector<MaterialName>>("tensor_materials")),
    _materials(getParam<std::vector<MaterialName>>("materials")),
    _def_grad_name(getParam<MaterialName>("def_grad_name")),
    _recover_num(getParam<unsigned int>("recover_num"))
{
}

// Function to add in auxvariables
void
RecoverVariablesAction::addAuxVariable(std::string name, InputParameters var_params, bool tensor)
{
  std::vector<std::string> conv = {"x", "y", "z"};

  if (tensor)
    for (unsigned int qp = 1; qp <= 8; qp++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
        {
          _problem->addAuxVariable("MooseVariable",
                                   name + "_" + conv[j] + conv[k] + "_" + std::to_string(qp),
                                   var_params);
        }
  if (!tensor)
    // Assuming 8 QPs
    // Starting at 1 because qps start at one in Sierra
    for (int qp = 1; qp <= 8; qp++)
      _problem->addAuxVariable("MooseVariable", name + "_" + std::to_string(qp), var_params);
}

void
RecoverVariablesAction::addAuxKernel(std::string var_name, std::string mat_name, bool tensor)
{
  std::vector<std::string> conv = {"x", "y", "z"};

  if (tensor)
    for (int qp = 1; qp <= 8; qp++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
        {
          InputParameters params = _factory.getValidParams("ADRankTwoAux");
          params.set<AuxVariableName>("variable") =
              var_name + "_" + conv[j] + conv[k] + "_" + std::to_string(qp);
          params.set<MaterialPropertyName>("rank_two_tensor") = mat_name;
          params.set<unsigned int>("selected_qp") = qp - 1;
          params.set<unsigned int>("index_i") = j;
          params.set<unsigned int>("index_j") = k;
          _problem->addAuxKernel("ADRankTwoAux",
                                 var_name + "_" + conv[j] + conv[k] + "_" + std::to_string(qp),
                                 params);
        }
  if (!tensor)
    for (int qp = 1; qp <= 8; qp++)
    {
      InputParameters params = _factory.getValidParams("ADMaterialRealAux");
      params.set<AuxVariableName>("variable") = var_name + "_" + std::to_string(qp);
      params.set<MaterialPropertyName>("property") = mat_name;
      params.set<unsigned int>("selected_qp") = qp - 1;
      _problem->addAuxKernel("ADMaterialRealAux", var_name + "_" + std::to_string(qp), params);
    }
}

void
RecoverVariablesAction::act()
{
  std::vector<std::string> conv = {"x", "y", "z"};

  if (_current_task == "add_aux_variable")
  {
    auto var_params = _factory.getValidParams("MooseVariable");
    var_params.set<MooseEnum>("order") = "CONSTANT";
    var_params.set<MooseEnum>("family") = "MONOMIAL";

    // Non tensor mats
    for (unsigned int i = 0; i < _materials.size(); i++)
      addAuxVariable(_materials[i], var_params, false);

    // Tensor mats
    for (unsigned int i = 0; i < _tensor_materials.size(); i++)
    {
      std::string matname = _tensor_materials[i];
      // Case for deformation gradients
      if (_tensor_materials[i].size() > 26)
        matname.erase(26, _tensor_materials[i].size() - 1);

      addAuxVariable(matname, var_params, true);
    }

    // Def grad
    for (unsigned int i = 0; i <= _recover_num; i++)
      addAuxVariable(_def_grad_name + "_" + std::to_string(i), var_params, true);
  }

  if (_current_task == "add_aux_kernel")
  {
    // Non tensor
    for (unsigned int i = 0; i < _materials.size(); i++)
      addAuxKernel(_materials[i], _materials[i], false);

    // tensor
    for (unsigned int i = 0; i < _tensor_materials.size(); i++)
    {
      std::string matname = _tensor_materials[i];

      if (_tensor_materials[i].size() > 26)
        matname.erase(26, _tensor_materials[i].size() - 1);
      addAuxKernel(matname, _tensor_materials[i], true);
    }

    std::vector<std::string> conv = {"x", "y", "z"};
    // def grad
    for (unsigned int i = 0; i < _recover_num; i++)
      addAuxKernel(
          _def_grad_name + "_" + std::to_string(i), _def_grad_name + "_" + std::to_string(i), true);

    addAuxKernel(_def_grad_name + "_" + std::to_string(_recover_num),
                 _def_grad_name + "_" + std::to_string(_recover_num),
                 true);
  }
}
