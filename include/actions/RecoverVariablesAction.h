#include "Action.h"
#include "InputParameters.h"
#include "MooseTypes.h"
#include "MooseEnum.h"
#include <unordered_map>

class RecoverVariablesAction : public Action

{
public:
  static InputParameters validParams();

  RecoverVariablesAction(const InputParameters & params);

  virtual void act() override;

protected:
  std::vector<MaterialName> _tensor_materials;

  std::vector<MaterialName> _materials;

  std::string _output_name;

  Real _qpnum;

  /// Element types
  const enum class Element { TET4_2nd, TET4_4th, TET10_4th, HEX8_3rd } _element;

private:
  const std::unordered_map<int, int> * _lookup;
};
