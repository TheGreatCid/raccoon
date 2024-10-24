#include "Action.h"
#include "InputParameters.h"
#include "MooseTypes.h"
#include "MooseEnum.h"

class RecoverVariablesAction : public Action

{
public:
  static InputParameters validParams();

  RecoverVariablesAction(const InputParameters & params);

  virtual void act() override;

protected:
  void addAuxVariable(std::string name, InputParameters var_params, bool tensor);

  void addAuxKernel(std::string var_name, std::string mat_name, bool tensor);

  std::vector<MaterialName> _tensor_materials;

  std::vector<MaterialName> _materials;

  MaterialName _def_grad_name;

  unsigned int _recover_num;
};
