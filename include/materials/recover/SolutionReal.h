//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "BaseNameInterface.h"

class SolutionUserObjectQP;

/**
 * A scalar material property read at each quadrature point from a recovery file
 */
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
  /// Not read. Requesting the old value makes the property stateful, which is what gets
  /// initStatefulProperties() called and keeps the recovered value across time steps.
  const MaterialProperty<Real> & _mat_old;

  const SolutionUserObjectQP * _solution_object_ptr;
};
