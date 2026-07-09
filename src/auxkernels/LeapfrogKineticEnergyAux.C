//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "LeapfrogKineticEnergyAux.h"

registerMooseObject("raccoonApp", LeapfrogKineticEnergyAux);

InputParameters
LeapfrogKineticEnergyAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Kinetic energy density in the leapfrog product form 1/2 rho (v^{n+1/2} . v^{n-1/2}), the "
      "quantity conserved by the central-difference integrator.  Couple the half-step velocity "
      "AuxVariables (v^{n+1/2}); their OLD state supplies v^{n-1/2}.  Density is read from a "
      "material property so it may vary in space.");
  params.addRequiredCoupledVar("velocity_x",
                               "X component of the leapfrog half-step velocity (v^{n+1/2}).");
  params.addRequiredCoupledVar("velocity_y",
                               "Y component of the leapfrog half-step velocity (v^{n+1/2}).");
  params.addRequiredCoupledVar("velocity_z",
                               "Z component of the leapfrog half-step velocity (v^{n+1/2}).");
  params.addParam<MaterialPropertyName>(
      "density", "density", "Name of the material property defining the density rho.");
  params.addParam<std::string>("base_name", "Mechanical property base name");
  return params;
}

LeapfrogKineticEnergyAux::LeapfrogKineticEnergyAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _density(getMaterialProperty<Real>(_base_name + "density")),
    _vel_x(coupledValue("velocity_x")),
    _vel_y(coupledValue("velocity_y")),
    _vel_z(coupledValue("velocity_z")),
    _vel_x_old(coupledValueOld("velocity_x")),
    _vel_y_old(coupledValueOld("velocity_y")),
    _vel_z_old(coupledValueOld("velocity_z"))
{
}

Real
LeapfrogKineticEnergyAux::computeValue()
{
  return 0.5 * _density[_qp] *
         (_vel_x[_qp] * _vel_x_old[_qp] + _vel_y[_qp] * _vel_y_old[_qp] +
          _vel_z[_qp] * _vel_z_old[_qp]);
}
