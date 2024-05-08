#pragma once

#include "ADTimeKernel.h"
#include "Material.h"
#include "SolutionUserObject.h"

// Forward Declarations
class TimeIntegrator;

// Parent class directly set to ADTimeKernel
class ADInertialForceRecover : public ADTimeKernel
{
public:
  static InputParameters validParams();

  ADInertialForceRecover(const InputParameters & parameters);

protected:
  virtual GenericReal<true> computeQpResidual();
  virtual Real computeQpJacobian();
  virtual void computeResidualAdditional();

private:
  const GenericMaterialProperty<Real, true> & _density;
  const VariableValue * _u_old;
  const VariableValue * _vel_old;
  const VariableValue * _accel_old;
  const bool _has_beta;
  const bool _has_gamma;
  const Real _beta;
  const Real _gamma;
  const bool _has_velocity;
  const bool _has_acceleration;
  const GenericMaterialProperty<Real, true> & _eta;
  const MaterialProperty<Real> & _density_scaling;
  const Real _alpha;

  // Velocity and acceleration calculated by time integrator
  const VariableValue * _u_dot_factor_dof;
  const VariableValue * _u_dotdot_factor_dof;
  const VariableValue * _u_dot_factor;
  const VariableValue * _u_dotdot_factor;
  const VariableValue * _u_dot_old;
  const VariableValue * _du_dot_du;
  const VariableValue * _du_dotdot_du;

  /// The TimeIntegrator
  TimeIntegrator & _time_integrator;

  const SolutionUserObject * _solution_object_ptr;
  VariableName _inert_name;
  VariableName _vel_old_name;
  VariableName _accel_old_name;

  // Use the variables from ADTimeKernel
  using ADTimeKernel::_dt;
  using ADTimeKernel::_i;
  using ADTimeKernel::_phi;
  using ADTimeKernel::_qp;
  using ADTimeKernel::_sys;
  using ADTimeKernel::_test;
  using ADTimeKernel::_u;
  using ADTimeKernel::_var;
};
