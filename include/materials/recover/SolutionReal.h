//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "InputParameters.h"
#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"
#include "MaterialProperty.h"
#include "MooseArray.h"
#include "SolutionUserObject.h"
#include "Qp_Mapping.h"

/**
 * Per-QP scalar material that reads values from a SolutionUserObject (typical
 * use: the recovery dump).  Templated on `is_ad` so the declared property can
 * be either AD (`SolutionReal`) for AD-consumer materials, or non-AD
 * (`SolutionRealNonAD`) for non-AD consumers like the framework MassMatrix
 * kernel that errors when a property has been declared as the wrong AD-ness.
 */
template <bool is_ad>
class SolutionRealTempl : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  SolutionRealTempl(const InputParameters & parameters);

  void initStatefulProperties(unsigned int n_points) override;

  void initialSetup() override;

protected:
  const std::string _mat_name;

  GenericMaterialProperty<Real, is_ad> & _mat;
  const MaterialProperty<Real> & _mat_old;

  const SolutionUserObject * _solution_object_ptr;

  QpMapping::Element _element;

  unsigned int _qpnum;

private:
  const std::unordered_map<int, int> * _lookup;
};

typedef SolutionRealTempl<true> SolutionReal;
typedef SolutionRealTempl<false> SolutionRealNonAD;
