//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ExplicitMixedOrder.h"

class SolutionUserObject;

/**
 * Recovery-aware variant of ExplicitMixedOrder (central-difference explicit
 * time integrator).
 *
 * The stock ExplicitMixedOrder::computeICs() initializes the integrator's
 * velocity vector (_sys.solutionUDot()) to zero at INITIAL via the forward
 * Euler estimate (solution - solution_old) / dt.  That's correct for a fresh
 * simulation start, but breaks the recover/restart pattern where the restart
 * needs to pick up at the non-zero velocity / acceleration the reference had
 * at dump_time.
 *
 * This class:
 *   - looks up a SolutionUserObject pointing at the reference's exodus dump,
 *   - after the base class init() finishes (so the default zero-IC has run),
 *     iterates local nodes and writes the recovered nodal vel / accel values
 *     directly into the integrator's solutionUDot / solutionUDotDot vectors,
 *   - leaves everything else (mass matrix, RHS evaluation, central-difference
 *     solve loop) unchanged.
 *
 * Input usage:
 *   [TimeIntegrator]
 *     type = ExplicitMixedOrderRecover
 *     mass_matrix_tag = 'mass'
 *     second_order_vars = 'disp_x disp_y disp_z'
 *     solution   = epsol                         # SolutionUserObject of dump
 *     disp_vars  = 'disp_x disp_y disp_z'        # nonlinear vars on restart side
 *     vel_vars   = 'vel_x vel_y vel_z'           # nodal aux var names in the dump
 *     accel_vars = 'accel_x accel_y accel_z'
 *   []
 *
 * The disp_vars / vel_vars / accel_vars lists are parallel: index i of each
 * list specifies the disp variable on the restart side and the corresponding
 * nodal vel / accel AuxVariable in the dump's exodus.
 *
 * Notes:
 *   - The restart's nonlinear disp variables themselves stay at zero at INITIAL
 *     (incremental displacement on the dump-deformed mesh).  This class only
 *     seeds the time-derivative vectors.
 *   - The dump's exodus (read by the SolutionUserObject) and the restart's
 *     mesh must be the same mesh from MOOSE's perspective.  In the typical
 *     remesh workflow, an external mapper has already projected the
 *     reference-side fields onto the restart's new mesh and written them to
 *     an exodus file whose nodal IDs match the restart's mesh; MOOSE does no
 *     interpolation on its side.  directValue(node, var_name) is then a
 *     direct nodal lookup.
 */
class ExplicitMixedOrderRecover : public ExplicitMixedOrder
{
public:
  static InputParameters validParams();

  ExplicitMixedOrderRecover(const InputParameters & parameters);

  virtual void init() override;
  virtual void preSolve() override;

protected:
  /// Iterate local nodes of each disp variable and populate solutionUDot /
  /// solutionUDotDot from the dump.  Called once at INITIAL after the base
  /// class init.
  void recoverVelAccel();

  /// The SolutionUserObject reading the reference dump.  Bound in init()
  /// (not in the constructor) because the UO has to exist by then.
  const SolutionUserObject * _sol_uo;

  /// Parallel lists: index i specifies the disp variable on the restart side
  /// (target DOFs) and the nodal vel / accel AuxVariable names in the dump.
  const std::vector<VariableName> _disp_var_names;
  const std::vector<std::string> _vel_var_names;
  const std::vector<std::string> _accel_var_names;

  /// Flips to false after the first preSolve().  Used to force
  /// _dt_old = _dt on the very first restart step so the central-difference
  /// velocity update (which scales accel by (_dt + _dt_old)/2) uses the
  /// full accel contribution -- matching what the reference's continued
  /// run would compute at the equivalent step, where _dt_old = previous dt.
  bool _first_solve_pending;
};
