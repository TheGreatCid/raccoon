//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "LumpedKineticEnergy.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"
#include "libmesh/numeric_vector.h"

registerMooseObject("raccoonApp", LumpedKineticEnergy);

InputParameters
LumpedKineticEnergy::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addClassDescription(
      "Kinetic energy in the integrator's own lumped-mass norm, 1/2 sum_i m_i (u_dot_i)^2, using "
      "the ExplicitMixedOrder 'mass_matrix_lumped' vector and solutionUDot.  This is the KE the "
      "central-difference scheme conserves against; comparing it between a reference and a recover "
      "run isolates the lumped-mass reconstruction error.");
  params.addParam<std::string>(
      "mass_vector", "mass_matrix_lumped", "Name of the integrator's lumped mass diagonal vector.");
  return params;
}

LumpedKineticEnergy::LumpedKineticEnergy(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _nl(_fe_problem.getNonlinearSystemBase(0)),
    _mass_vector_name(getParam<std::string>("mass_vector")),
    _ke(0.0)
{
}

void
LumpedKineticEnergy::execute()
{
  if (!_nl.hasVector(_mass_vector_name))
    mooseError("LumpedKineticEnergy: the nonlinear system has no vector named '",
               _mass_vector_name,
               "'.  This postprocessor requires an ExplicitMixedOrder-family time integrator.");

  const NumericVector<Number> & m = _nl.getVector(_mass_vector_name);

  const NumericVector<Number> * const v_ptr = _nl.solutionUDot();
  if (!v_ptr)
    mooseError("LumpedKineticEnergy: solutionUDot is not available.");
  const NumericVector<Number> & v = *v_ptr;

  // KE = 1/2 sum_i m_i v_i^2 over locally-owned dofs, then reduce.
  Real ke = 0.0;
  for (dof_id_type i = v.first_local_index(); i < v.last_local_index(); ++i)
  {
    const Real vi = v(i);
    ke += 0.5 * m(i) * vi * vi;
  }
  _communicator.sum(ke);
  _ke = ke;
}

PostprocessorValue
LumpedKineticEnergy::getValue() const
{
  return _ke;
}
