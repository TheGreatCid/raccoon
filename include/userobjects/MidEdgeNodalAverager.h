//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "GeneralUserObject.h"

/**
 * For each second-order nodal AuxVariable in `variables`, replace the value at
 * every mid-edge node with the arithmetic mean of the two corner nodes of that
 * edge.  Restart-time fix for cross-mesh recovery: SolutionIC samples the old
 * mesh pointwise at new-mesh node positions, which can leave mid-edge values
 * inconsistent with their corner neighbors and excite the inertia residual via
 * the (1/(beta*dt)) * v_old term.
 */
class MidEdgeNodalAverager : public GeneralUserObject
{
public:
  static InputParameters validParams();
  MidEdgeNodalAverager(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual void finalize() override {}

protected:
  std::vector<AuxVariableName> _var_names;
};
