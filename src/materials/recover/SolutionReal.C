#include "BaseNameInterface.h"
#include "InputParameters.h"
#include "Material.h"
#include "MooseTypes.h"
#include "RankTwoTensorForward.h"
#include "SolutionUserObject.h"
#include "SolutionReal.h"
#include "Qp_Mapping.h"

registerADMooseObject("raccoonApp", SolutionReal);

InputParameters
SolutionReal::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();

  params.addClassDescription("Getting tensor from solution object");
  params.addParam<UserObjectName>("solution", "The SolutionUserObject to extract data from.");
  params.addParam<Real>("num_qps", 8, "Number of QPs");
  params.addParam<std::string>("mat_name", "stress", "Name of tensor");
  params.addRequiredParam<MooseEnum>(
      "element", MooseEnum("TET4_2nd TET4_4th TET10_4th HEX8_3rd"), "The element type");
  return params;
}

SolutionReal::SolutionReal(const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    _mat_name(getParam<std::string>("mat_name")),
    _mat(declareADProperty<Real>(prependBaseName(_mat_name + "_sol"))),
    _mat_old(getMaterialPropertyOld<Real>(prependBaseName(_mat_name + "_sol"))),
    _solution_object_ptr(NULL),
    _qpnum(getParam<Real>("num_qps")),
    _element(getParam<MooseEnum>("element").getEnum<Element>())
{
  // _mesh->elemTypes();
  // Picking mapping for MOOSE to SIERRA qp numbering
  switch (_element)
  {
    // _mesh->elemTypes();
    // Picking mapping for MOOSE to SIERRA qp numbering
    case Element::TET4_2nd:
      _lookup = &Qp_Mapping::TET4_2nd_lookup_rev;
      _qpnum = 4;
      break;
    case Element::TET4_4th:
      _lookup = &Qp_Mapping::TET4_4th_lookup_rev;
      _qpnum = 5;
      break;
    case Element::TET10_4th:
      _lookup = &Qp_Mapping::TET10_4th_lookup_rev;
      _qpnum = 11;
      break;
    case Element::HEX8_3rd:
      _lookup = &Qp_Mapping::HEX8_3rd_lookup_rev;
      _qpnum = 8;
      break;
    default:
      mooseError("Unknown element type");
  }
}

void
SolutionReal::initialSetup()
{
  _solution_object_ptr = &getUserObject<SolutionUserObject>("solution");
}

void
SolutionReal::initStatefulProperties(unsigned int n_points)
{
  unsigned int qp_max = _qpnum;

  auto formatQP = [qp_max](unsigned int qp)
  {
    if (qp_max < 10)
      return std::to_string(qp); // Single digit
    else
      return (qp < 10) ? "0" + std::to_string(qp) : std::to_string(qp); // Two digits
  };

  for (_qp = 0; _qp < n_points; ++_qp)
  {
    Real _qp_sel = _lookup->find(_qp + 1)->second;
    // Populate from solution user object
    _mat[_qp] = _solution_object_ptr->pointValue(
        _t, _current_elem->true_centroid(), _mat_name + "_" + formatQP(_qp_sel), nullptr);
  }
}