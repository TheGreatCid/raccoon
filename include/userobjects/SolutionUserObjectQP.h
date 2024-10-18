#pragma once

#include "SolutionUserObject.h"

class SolutionUserObjectQP : public SolutionUserObject
{
public:
  static InputParameters validParams();
  SolutionUserObjectQP(const InputParameters & parameters);

protected:
  std::vector<MaterialName> _tensor_materials;
  std::vector<MaterialName> _materials;
};
