//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "LargeDeformationPlasticityModel.h"

/// Intermediate base for J2 plasticity models whose state can be recovered at the quadrature
/// points. Its initQpStatefulProperties seeds _ep and _Fp from the solution when recover = true.
class LargeDeformationJ2PlasticityBase : public LargeDeformationPlasticityModel
{
public:
  static InputParameters validParams();

  LargeDeformationJ2PlasticityBase(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;

  /// A non-negative scalar property of the recovered solution at the current quadrature point
  /// (remapping can leave small negative values)
  Real recoverNonNegative(const std::string & property) const;

  /// A tensor property of the recovered solution at the current quadrature point
  RankTwoTensor recoverTensor(const std::string & property) const;
};
