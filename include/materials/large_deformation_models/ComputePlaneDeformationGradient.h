//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ComputeDeformationGradient.h"

/**
 * Plane-stress deformation gradient. Identical to ComputeDeformationGradient (including the full
 * recover / F-bar / polar-decomposition machinery) except that the out-of-plane stretch is supplied
 * by a coupled out-of-plane strain variable: F_zz = 1 + out_of_plane_strain, imposed at every
 * quadrature point rather than only for the axisymmetric (RZ) case.
 */
class ComputePlaneDeformationGradient : public ComputeDeformationGradient
{
public:
  static InputParameters validParams();

  ComputePlaneDeformationGradient(const InputParameters & parameters);

protected:
  ADReal computeQpOutOfPlaneGradDisp() override;
  void applyOutOfPlaneGradDisp(ADRankTwoTensor & A) override;

  /// Out-of-plane strain variable driving F_zz = 1 + out_of_plane_strain (plane stress).
  const ADVariableValue & _out_of_plane_strain;
};
