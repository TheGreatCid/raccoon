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
  params.addParam<std::string>("output_name", "exodusqp", "name of output for variables");
  return params;
}

RecoverVariablesAction::RecoverVariablesAction(const InputParameters & params)
  : Action(params),
    _tensor_materials(getParam<std::vector<MaterialName>>("tensor_materials")),
    _materials(getParam<std::vector<MaterialName>>("materials")),
    _output_name(getParam<std::string>("output_name"))
{
}

void
RecoverVariablesAction::act()
{
  auto dim = _mesh->dimension();
  unsigned int qp_max = 2;

  std::vector<std::string> conv = {"x", "y", "z"};

  if (_current_task == "add_aux_variable")
  {
    auto var_params = _factory.getValidParams("MooseVariable");
    var_params.set<MooseEnum>("order") = "CONSTANT";
    var_params.set<MooseEnum>("family") = "MONOMIAL";
    var_params.set<std::vector<OutputName>>("outputs") = {_output_name};
    // Non tensor mats
    for (unsigned int i = 0; i < _materials.size(); i++)
    {
      // Assuming 8 QPs
      // Starting at 1 because qps start at one in Sierra
      for (unsigned int qp = 1; qp <= qp_max; qp++)
        _problem->addAuxVariable(
            "MooseVariable", _materials[i] + "_" + std::to_string(qp), var_params);
    }

    // Tensor mats
    for (unsigned int i = 0; i < _tensor_materials.size(); i++)
    {
      // Assuming 8 QPs
      // Starting at 1 because qps start at one in Sierra
      for (unsigned int qp = 1; qp <= qp_max; qp++)
        for (unsigned int j = 0; j < dim; j++)
          for (unsigned int k = 0; k < dim; k++)
          {
            std::string matname = _tensor_materials[i];
            if (_tensor_materials[i].size() > 26)
              matname.erase(26, _tensor_materials[i].size() - 1);
            _problem->addAuxVariable("MooseVariable",
                                     matname + "_" + conv[j] + conv[k] + "_" + std::to_string(qp),
                                     var_params);
          }
    }
  }

  if (_current_task == "add_aux_kernel")
  {
    // Non tensor
    for (unsigned int i = 0; i < _materials.size(); i++)
      // assuming 8 QPs
      for (unsigned int qp = 1; qp <= qp_max; qp++)
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
    {
      std::string matname = _tensor_materials[i];
      if (_tensor_materials[i].size() > 26)
        matname.erase(26, _tensor_materials[i].size() - 1);
      // assuming 8 QPs
      for (unsigned int qp = 1; qp <= qp_max; qp++)
        for (unsigned int j = 0; j < dim; j++)
          for (unsigned int k = 0; k < dim; k++)
          {
            InputParameters params = _factory.getValidParams("ADRankTwoAux");
            params.set<AuxVariableName>("variable") =
                matname + "_" + conv[j] + conv[k] + "_" + std::to_string(qp);
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
}
