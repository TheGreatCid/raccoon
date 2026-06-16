#pragma once

#include "DirichletBCBase.h"
#include "Function.h"

/**
 * Spatially-aware version of MOOSE's PresetDisplacement.
 *
 * The stock PresetDisplacement (modules/solid_mechanics/src/bcs/PresetDisplacement.C)
 * calls `_function.timeDerivative(_t)` -- the helper at
 * `framework/include/functions/Function.h:220-223` that forwards to
 * `timeDerivative(t, p)` with a default-constructed `p = (0,0,0)`.  For any
 * spatially-varying BC the spatial factor vanishes at the origin and the BC
 * contributes nothing.
 *
 * This variant evaluates the function's time derivative at the
 *   - current BC node's position when no `coupled_x` (etc.) AuxVariable is
 *     supplied (typical reference-run use case where the node's coordinates
 *     ARE the undeformed coordinates), or
 *   - the value of the supplied `coupled_x`/`coupled_y`/`coupled_z`
 *     AuxVariables at the node (typical restart use case where the node sits
 *     on the dump-deformed mesh and the original undeformed coordinates are
 *     carried as recovered AuxVariables like X0).
 *
 * Everything else (the Newmark forward-update formula for the enforced
 * displacement) is identical to PresetDisplacement, so the BC is Newmark-
 * consistent and does not suffer the 1/(beta*dt) accel-blowup of an
 * ADMatchedValueBC + NewmarkAccelAux back-computation when the vel/accel
 * AuxVariables hold values that don't perfectly match the function's
 * analytical derivatives.
 */
class PresetDisplacementSpatial : public DirichletBCBase
{
public:
  static InputParameters validParams();

  PresetDisplacementSpatial(const InputParameters & parameters);

protected:
  virtual Real computeQpValue() override;

  /// Build the Point at which to evaluate the driving function:
  /// AuxVariable-supplied coord for any axis that's been coupled, falling
  /// back to the node's own coordinate for the others.
  Point spatialPoint() const;

  const VariableValue & _u_old;
  const Real _scale_factor;
  const Function & _function;
  const VariableValue & _vel_old;
  const VariableValue & _accel_old;
  const Real _beta;

  const bool _has_coupled_x;
  const bool _has_coupled_y;
  const bool _has_coupled_z;
  const VariableValue * const _coupled_x;
  const VariableValue * const _coupled_y;
  const VariableValue * const _coupled_z;
};
