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

class SolutionReal : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  SolutionReal(const InputParameters & parameters);

  void initStatefulProperties(unsigned int n_points) override;

  void initialSetup() override;

protected:
  const std::string _mat_name;

  ADMaterialProperty<Real> & _mat;
  const MaterialProperty<Real> & _mat_old;

  const SolutionUserObject * _solution_object_ptr;
  Real _qpnum;

  /// Element types
  const enum class Element { TET4_2nd, TET4_4th, TET10_4th, HEX8_3rd } _element;

private:
  const std::unordered_map<int, int> * _lookup;
};