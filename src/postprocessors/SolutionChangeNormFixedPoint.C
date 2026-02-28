//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "SolutionChangeNormFixedPoint.h"
#include "Moose.h"

registerMooseObject("raccoonApp", SolutionChangeNormFixedPoint);

InputParameters
SolutionChangeNormFixedPoint::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Computes the L2 norm of the change in a variable between consecutive fixed-point "
      "(Picard) iterations: ||d^(k+1) - d^(k)||_2. Intended for use as a staggered-scheme "
      "convergence criterion via PostprocessorConvergence.");
  params.addRequiredCoupledVar("variable", "The variable to monitor for fixed-point convergence");
  // Default to running at the fixed-point convergence check, after each sub-app solve and
  // damage transfer, but before the convergence decision is made.
  params.set<ExecFlagEnum>("execute_on") = {EXEC_MULTIAPP_FIXED_POINT_CONVERGENCE};
  return params;
}

SolutionChangeNormFixedPoint::SolutionChangeNormFixedPoint(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _fp_var(nullptr)
{
  auto * base_var = getVar("variable", 0);
  _fp_var = dynamic_cast<MooseVariableFE<Real> *>(base_var);
  if (!_fp_var)
    mooseError(name(),
               ": variable '",
               base_var->name(),
               "' must be a standard (Lagrange) finite element variable.");
}

void
SolutionChangeNormFixedPoint::initialize()
{
  ElementIntegralPostprocessor::initialize();

  auto & sys = _fp_var->sys();

  // solutionState(0) is the current solution, updated by the MultiAppCopyTransfer to d^(k).
  // We must localize it directly — coupledValue() reads current_local_solution which is NOT
  // synced after a transfer writes to solution(), so it would still contain d^(k-1).
  sys.solutionState(0).localize(_current_local);

  // solutionState(1, FixedPoint) is the PREVIOUS_FP_SOLUTION_TAG vector — d^(k-1), copied
  // there by SystemBase::copyPreviousFixedPointSolutions() before the current solve began.
  // MOOSE only calls needSolutionState on _solver_sys, not _aux, so hasSolutionState always
  // returns false here. Call the non-const overload directly to lazily create the state.
  // On the first FP iteration ever the state initialises to zero; subsequent iterations are
  // populated correctly by copyPreviousFixedPointSolutions() once the state exists.
  sys.solutionState(1, Moose::SolutionIterationType::FixedPoint).localize(_fp_old_local);
}

Real
SolutionChangeNormFixedPoint::computeQpIntegral()
{
  // Interpolate the previous fixed-point solution at the current quadrature point.
  // dofIndices() and phi() are populated for the current element during reinit().
  const auto & dof_indices = _fp_var->dofIndices();
  const auto & phi = _fp_var->phi();

  Real u_curr = 0.0, u_fp_old = 0.0;
  for (unsigned int i = 0; i < dof_indices.size(); ++i)
  {
    u_curr   += _current_local[dof_indices[i]] * phi[i][_qp];
    u_fp_old += _fp_old_local[dof_indices[i]] * phi[i][_qp];
  }

  const Real diff = u_curr - u_fp_old;
  return diff * diff;
}

void
SolutionChangeNormFixedPoint::finalize()
{
  ElementIntegralPostprocessor::finalize();
  if (processor_id() == 0)
    Moose::out << "  " << name() << " = " << std::sqrt(_integral_value) << "\n";
}

Real
SolutionChangeNormFixedPoint::getValue() const
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}
