#include "Factory.h"
#include "FEProblem.h"
#include "AddVariableAction.h"
#include "InputParameters.h"
#include "MooseUtils.h"
#include "NonlinearSystemBase.h"
#include "Parser.h"
#include "Qp_Mapping.h"
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
  params.addRequiredParam<MooseEnum>(
      "element", MooseEnum(QpMapping::ELEMENT_ENUM_DEFINITION), "The element type");
  params.addParam<std::vector<SubdomainName>>(
      "block",
      "The list of blocks (ids or names) the generated AuxKernels are restricted to. If "
      "specified, the AuxKernel names are suffixed with the block names.");

  // Defaults to the stock AuxKernel schedule so existing inputs are unchanged.  This block
  // generates one AuxKernel per tensor component per QP (e.g. 4 tensors x 9 components x 11 QPs
  // = 396 kernels for TET10_4th), and on the default schedule every one of them re-evaluates its
  // material property on each linear iteration.  Since the recovery variables are only ever
  // consumed by the dump file, narrowing this to the schedule the output actually uses (e.g.
  // 'FINAL', or the sync times the recovery file is written at) removes that cost entirely.
  ExecFlagEnum exec_enum = MooseUtils::getDefaultExecFlagEnum();
  exec_enum.addAvailableFlags(EXEC_PRE_DISPLACE);
  exec_enum = {EXEC_LINEAR, EXEC_TIMESTEP_END};
  params.addParam<ExecFlagEnum>(
      "execute_on", exec_enum, "When the generated recovery AuxKernels are evaluated.");
  return params;
}

RecoverVariablesAction::RecoverVariablesAction(const InputParameters & params)
  : Action(params),
    _tensor_materials(getParam<std::vector<MaterialName>>("tensor_materials")),
    _materials(getParam<std::vector<MaterialName>>("materials")),
    _output_name(getParam<std::string>("output_name")),
    _element(getParam<MooseEnum>("element").getEnum<QpMapping::Element>()),
    _blocks(isParamValid("block") ? getParam<std::vector<SubdomainName>>("block")
                                  : std::vector<SubdomainName>()),
    _execute_on(getParam<ExecFlagEnum>("execute_on"))
{
  _lookup = QpMapping::getLookup(_element, _qpnum, /*reversed=*/false);
}

void
RecoverVariablesAction::act()
{
  // Tensors are always dumped with all 3x3 components, which is what the recovery readers expect
  const unsigned int dim = 3;

  // The variables are named after the quadrature points of the file, numbered from 1; _lookup
  // gives the MOOSE quadrature point each of them holds.
  if (_current_task == "add_aux_variable")
  {
    auto var_params = _factory.getValidParams("MooseVariable");
    var_params.set<MooseEnum>("order") = "CONSTANT";
    var_params.set<MooseEnum>("family") = "MONOMIAL";
    var_params.set<std::vector<OutputName>>("outputs") = {_output_name};

    for (unsigned int qp = 1; qp <= _qpnum; qp++)
    {
      for (const auto & material : _materials)
        _problem->addAuxVariable(
            "MooseVariable", QpMapping::variableName(material, qp, _qpnum), var_params);
      for (const auto & material : _tensor_materials)
        for (unsigned int j = 0; j < dim; j++)
          for (unsigned int k = 0; k < dim; k++)
            _problem->addAuxVariable(
                "MooseVariable", QpMapping::variableName(material, j, k, qp, _qpnum), var_params);
    }
  }

  if (_current_task == "add_aux_kernel")
  {
    // Suffix appended to AuxKernel names when block restricted, e.g. "_block1_block2"
    std::string block_suffix;
    for (const auto & block : _blocks)
      block_suffix += "_" + block;

    for (unsigned int qp = 1; qp <= _qpnum; qp++)
    {
      const unsigned int moose_qp = QpMapping::getQP(qp, _lookup) - 1;

      for (const auto & material : _materials)
      {
        const auto variable = QpMapping::variableName(material, qp, _qpnum);
        InputParameters params = _factory.getValidParams("ADMaterialRealAux");
        params.set<AuxVariableName>("variable") = variable;
        params.set<MaterialPropertyName>("property") = material;
        params.set<unsigned int>("selected_qp") = moose_qp;
        params.set<ExecFlagEnum>("execute_on") = _execute_on;
        if (!_blocks.empty())
          params.set<std::vector<SubdomainName>>("block") = _blocks;
        _problem->addAuxKernel("ADMaterialRealAux", variable + block_suffix, params);
      }

      for (const auto & material : _tensor_materials)
        for (unsigned int j = 0; j < dim; j++)
          for (unsigned int k = 0; k < dim; k++)
          {
            const auto variable = QpMapping::variableName(material, j, k, qp, _qpnum);
            InputParameters params = _factory.getValidParams("ADRankTwoAux");
            params.set<AuxVariableName>("variable") = variable;
            params.set<MaterialPropertyName>("rank_two_tensor") = material;
            params.set<unsigned int>("selected_qp") = moose_qp;
            params.set<unsigned int>("index_i") = j;
            params.set<unsigned int>("index_j") = k;
            params.set<ExecFlagEnum>("execute_on") = _execute_on;
            if (!_blocks.empty())
              params.set<std::vector<SubdomainName>>("block") = _blocks;
            _problem->addAuxKernel("ADRankTwoAux", variable + block_suffix, params);
          }
    }
  }
}
