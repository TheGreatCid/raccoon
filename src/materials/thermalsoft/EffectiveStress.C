//* David Development

#include "EffectiveStress.h"
#include "RaccoonUtils.h"

registerMooseObject("raccoonApp", EffectiveStress);

InputParameters
EffectiveStress::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates effective stress");
  return params;
}

EffectiveStress::EffectiveStress(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _stress(getADMaterialPropertyByName<RankTwoTensor>("stress")),
    _effective_stress(declareADProperty<ADReal>(prependBaseName("effective_stress")))
{
}

void
EffectiveStress::initQpStatefulProperties()
{
  _effective_stress[_qp] = 0;
}

void
EffectiveStress::computeQpProperties()
{
  ADRankTwoTensor devstress = _stress[_qp].deviatoric();
  _effective_stress[_qp] = std::sqrt((3 / 2) * devstress.doubleContraction(devstress));
}
