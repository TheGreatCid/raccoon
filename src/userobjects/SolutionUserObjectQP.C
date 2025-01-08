#include "SolutionUserObjectQP.h"

registerMooseObject("raccoonApp", SolutionUserObjectQP);

InputParameters
SolutionUserObjectQP::validParams()
{
  InputParameters params = SolutionUserObject::validParams();
  params.addParam<std::vector<MaterialName>>("tensor_materials", "materials to output qps on");
  params.addParam<std::vector<MaterialName>>("materials", "materials to output qps on");
  return params;
}

SolutionUserObjectQP::SolutionUserObjectQP(const InputParameters & parameters)
  : SolutionUserObject(parameters),
    _tensor_materials(getParam<std::vector<MaterialName>>("tensor_materials")),
    _materials(getParam<std::vector<MaterialName>>("materials"))
{
  unsigned int qp_max = 8;
  int dim = 3;

  if (!_tensor_materials.empty())
  {
    std::vector<std::string> conv = {"x", "y", "z"};
    for (unsigned int i = 0; i < _tensor_materials.size(); i++)
    {
      std::string matname = _tensor_materials[i];
      if (_tensor_materials[i].size() > 26)
        matname.erase(26, _tensor_materials[i].size() - 1);
      for (unsigned int qp = 1; qp <= qp_max; qp++)
        for (int j = 0; j < dim; j++)
          for (int k = 0; k < dim; k++)
          {
            _system_variables.push_back(matname + "_" + conv[j] + conv[k] + "_" +
                                        std::to_string(qp));
          }
    }
  }
  if (!_materials.empty())
  {
    for (unsigned int i = 0; i < _materials.size(); i++)
    {
      // Assuming 8 QPs
      // Starting at 1 because qps start at one in Sierra
      for (unsigned int qp = 1; qp <= qp_max; qp++)
        _system_variables.push_back(_materials[i] + "_" + std::to_string(qp));
    }
  }
}
