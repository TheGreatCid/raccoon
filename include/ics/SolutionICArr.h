//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Material.h"

// #include "InitialCondition.h"

class SolutionUserObject;

/**
 * Class for reading an initial condition from a solution user object
 */
class SolutionICArr : public Material
{
public:
  SolutionICArr(const InputParameters & parameters);

  virtual void initialSetup() override;

  Real returnval(const int & j, const Point & p);

protected:
  /// SolutionUserObject containing the solution of interest
  const SolutionUserObject & _solution_object;

  /// The variable name extracted from the SolutionUserObject
  std::vector<std::string> _system_variables;

  /// Remapped IDs from the current mesh to the ExodusII mesh
  std::set<SubdomainID> _exo_block_ids;

  /// The system object
  SystemBase & _sys;

public:
  static InputParameters validParams();
};
