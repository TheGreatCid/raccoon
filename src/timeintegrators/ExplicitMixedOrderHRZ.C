//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ExplicitMixedOrderHRZ.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/nonlinear_solver.h"

registerMooseObject("raccoonApp", ExplicitMixedOrderHRZ);

InputParameters
ExplicitMixedOrderHRZ::validParams()
{
  InputParameters params = ExplicitMixedOrder::validParams();
  params.addClassDescription(
      "Variant of ExplicitMixedOrder that supports Hinton-Rock-Zienkiewicz (HRZ) "
      "mass-matrix lumping in addition to the parent's row-sum lumping.  HRZ is "
      "required for stability with TET10 elements where row-sum lumping produces "
      "near-zero or negative diagonal entries at corner nodes.");
  params.addParam<MooseEnum>(
      "mass_lumping_type",
      MooseEnum("row_sum hrz", "row_sum"),
      "Mass-matrix lumping method.  'row_sum' (default) reproduces the parent "
      "ExplicitMixedOrder behavior exactly: lumped_i = sum_j M_ij.  'hrz' uses "
      "the Hinton-Rock-Zienkiewicz scaling: take the diagonal entries of the "
      "consistent mass matrix (which are positive for any conforming basis) and "
      "rescale them so their sum equals the total mass.  HRZ is the standard "
      "fix for TET10 + central-difference where row-sum produces ill-conditioned "
      "lumped diagonals.");
  return params;
}

ExplicitMixedOrderHRZ::ExplicitMixedOrderHRZ(const InputParameters & parameters)
  : ExplicitMixedOrder(parameters),
    _mass_lumping_type(getParam<MooseEnum>("mass_lumping_type"))
{
}

void
ExplicitMixedOrderHRZ::solve()
{
  // Fast path: HRZ not requested -> exactly the parent's behavior.
  if (_mass_lumping_type == "row_sum")
  {
    ExplicitMixedOrder::solve();
    return;
  }

  // HRZ path.  Mirrors ExplicitMixedOrder::solve() but swaps the lumping step
  // for HRZ.  The lumping must happen BEFORE evaluateRHSResidual because
  // ExplicitFunctionDirichletBC's residual reads the lumped mass diagonal
  // directly -- if BC residual and integrator inverse-mass use different
  // lumping, the BC kinematics are inconsistent.
  auto mass_tag = massMatrixTagID();

  _n_nonlinear_iterations = 0;
  _n_linear_iterations = 0;

  _current_time = _fe_problem.time();

  auto & mass_matrix = _nonlinear_implicit_system->get_system_matrix();

  if (_mesh_changed)
    updateDOFIndices();

  if (_mesh_changed || !_constant_mass)
  {
    // Assemble the consistent mass matrix into the tagged matrix.
    _fe_problem.computeJacobianTag(
        *_nonlinear_implicit_system->current_local_solution, mass_matrix, mass_tag);

    // HRZ lumping:
    //   1) total mass = sum of all entries  (computed as sum of row-sums)
    //   2) overwrite the lumped vector with the matrix's actual diagonal
    //   3) scale diagonal so its sum equals total mass
    mass_matrix.vector_mult(*_mass_matrix_lumped, *_ones);
    const Real m_total = _mass_matrix_lumped->sum();

    mass_matrix.get_diagonal(*_mass_matrix_lumped);
    const Real diag_sum = _mass_matrix_lumped->sum();
    if (diag_sum <= 0.0)
      mooseError("ExplicitMixedOrderHRZ: HRZ lumping produced non-positive "
                 "diagonal sum (",
                 diag_sum,
                 ").  Mass matrix may have zero or negative diagonal entries.");
    _mass_matrix_lumped->scale(m_total / diag_sum);
    _mass_matrix_lumped->close();

    *_mass_matrix_diag_inverted = *_mass_matrix_lumped;
    _mass_matrix_diag_inverted->reciprocal();
    _mass_matrix_diag_inverted->close();
  }

  _mesh_changed = false;

  // Remainder is bit-identical to the parent: set time, update system,
  // evaluate residual, solve, write back.
  _fe_problem.time() = _fe_problem.timeOld();
  _nonlinear_implicit_system->update();

  evaluateRHSResidual();

  bool converged = performExplicitSolve(mass_matrix);
  _nl->overwriteNodeFace(*_nonlinear_implicit_system->solution);

  *_nonlinear_implicit_system->solution = _nl->solutionOld();
  *_nonlinear_implicit_system->solution += *_solution_update;

  _nonlinear_implicit_system->update();

  _nl->setSolution(*_nonlinear_implicit_system->current_local_solution);
  _nonlinear_implicit_system->nonlinear_solver->converged = converged;
}
