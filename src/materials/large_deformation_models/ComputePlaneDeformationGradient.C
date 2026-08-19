//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ComputePlaneDeformationGradient.h"

registerADMooseObject("raccoonApp", ComputePlaneDeformationGradient);

InputParameters
ComputePlaneDeformationGradient::validParams()
{
  InputParameters params = ComputeDeformationGradient::validParams();
  params.addClassDescription(
      "Computes the deformation gradient for a plane-stress problem: the out-of-plane stretch is "
      "F_zz = 1 + out_of_plane_strain, supplied by a coupled out-of-plane strain variable. Inherits "
      "the full recover / F-bar / polar-decomposition machinery of ComputeDeformationGradient.");
  params.addCoupledVar(
      "out_of_plane_strain",
      "The out-of-plane strain variable; sets F_zz = 1 + out_of_plane_strain. If not coupled it is "
      "zero, recovering plane strain (F_zz = 1).");
  return params;
}

ComputePlaneDeformationGradient::ComputePlaneDeformationGradient(const InputParameters & parameters)
  : ComputeDeformationGradient(parameters),
    _out_of_plane_strain(adCoupledValue("out_of_plane_strain"))
{
}

ADReal
ComputePlaneDeformationGradient::computeQpOutOfPlaneGradDisp()
{
  // Out-of-plane stretch F_zz = 1 + out_of_plane_strain (engineering out-of-plane strain).
  return 1.0 + _out_of_plane_strain[_qp];
}

void
ComputePlaneDeformationGradient::applyOutOfPlaneGradDisp(ADRankTwoTensor & A)
{
  // Plane stress: impose the out-of-plane stretch at every QP (not gated on RZ). A holds the
  // displacement-gradient tensor before the identity is added, so A(2,2) = F_zz - 1 = strain_zz.
  A(2, 2) = computeQpOutOfPlaneGradDisp() - 1.0;
}
