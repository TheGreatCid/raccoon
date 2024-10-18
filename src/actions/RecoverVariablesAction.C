#include "Factory.h"
#include "FEProblem.h"
#include "AddVariableAction.h"
#include "InputParameters.h"
#include "NonlinearSystemBase.h"
#include "Parser.h"
#include "RecoverVariablesAction.h"

registerMooseAction("raccoonApp", RecoverVariablesAction, "add_aux_variable");
registerMooseAction("raccoonApp", RecoverVariablesAction, "add_aux_kernel");

InputParameters
RecoverVariablesAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addParam<std::vector<MaterialName>>("tensor_materials", "materials to output qps on");
  params.addParam<std::vector<MaterialName>>("materials", "materials to output qps on");

  return params;
}

RecoverVariablesAction::RecoverVariablesAction(const InputParameters & params)
  : Action(params),
    _tensor_materials(getParam<std::vector<MaterialName>>("tensor_materials")),
    _materials(getParam<std::vector<MaterialName>>("materials"))
{
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
    {
      // Assuming 8 QPs
      // Starting at 1 because qps start at one in Sierra
      for (int qp = 1; qp <= 8; qp++)
        _problem->addAuxVariable(
            "MooseVariable", _materials[i] + "_" + std::to_string(qp), var_params);
    }

    // Tensor mats
    for (unsigned int i = 0; i < _tensor_materials.size(); i++)
    {
      // Assuming 8 QPs
      // Starting at 1 because qps start at one in Sierra
      for (unsigned int qp = 1; qp <= 8; qp++)
        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
          {
            _problem->addAuxVariable("MooseVariable",
                                     _tensor_materials[i] + "_" + conv[j] + conv[k] + "_" +
                                         std::to_string(qp),
                                     var_params);
            std::cout << _tensor_materials[i] + "_" + conv[j] + conv[k] + "_" + std::to_string(qp)
                      << std::endl;
          }
    }
  }

  if (_current_task == "add_aux_kernel")
  {
    // Non tensor
    for (unsigned int i = 0; i < _materials.size(); i++)
      // assuming 8 QPs
      for (int qp = 1; qp <= 8; qp++)
      {
        InputParameters params = _factory.getValidParams("ADMaterialRealAux");
        params.set<AuxVariableName>("variable") = _materials[i] + "_" + std::to_string(qp);
        params.set<MaterialPropertyName>("property") = _materials[i];
        params.set<unsigned int>("selected_qp") = qp - 1;
        _problem->addAuxKernel(
            "ADMaterialRealAux", _materials[i] + "_" + std::to_string(qp), params);
      }

    // tensor
    for (unsigned int i = 0; i < _tensor_materials.size(); i++)
      // assuming 8 QPs
      for (int qp = 1; qp <= 8; qp++)
        for (int j = 0; j < 3; j++)
          for (int k = 0; k < 3; k++)
          {
            InputParameters params = _factory.getValidParams("ADRankTwoAux");
            params.set<AuxVariableName>("variable") =
                _tensor_materials[i] + "_" + conv[j] + conv[k] + "_" + std::to_string(qp);
            params.set<MaterialPropertyName>("rank_two_tensor") = _tensor_materials[i];
            params.set<unsigned int>("selected_qp") = qp - 1;
            params.set<unsigned int>("index_i") = j;
            params.set<unsigned int>("index_j") = k;
            _problem->addAuxKernel("ADRankTwoAux",
                                   _tensor_materials[i] + "_" + conv[j] + conv[k] + "_" +
                                       std::to_string(qp),
                                   params);
          }
  }
}
