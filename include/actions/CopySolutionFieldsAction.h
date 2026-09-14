#pragma once

#include "Action.h"

/**
 * Copies the nodal and elemental variables of an Exodus file onto the current mesh as
 * AuxVariables, without listing them in the input file.
 *
 * The variable and block names are read from the file header, and for every variable the action
 * adds an AuxVariable (nodal: LAGRANGE of the mesh's order; elemental: CONSTANT MONOMIAL) and a
 * SolutionIC reading one time step through a shared SolutionUserObject. The values are looked up
 * by position, so the mesh must be in the file's configuration; it may have extra nodes (e.g. a
 * crack opened by BreakMeshByBlockGenerator) and different subdomains.
 */
class CopySolutionFieldsAction : public Action
{
public:
  static InputParameters validParams();

  CopySolutionFieldsAction(const InputParameters & params);

  virtual void act() override;

protected:
  /// Reads the variable and block names from the file header (once) and applies the filters
  void readHeader();

  /// Name of the AuxVariable holding a file variable
  std::string auxVariableName(const std::string & file_variable) const;

  /// The Exodus file to copy from
  const FileName & _file;

  /// Nodal variables to copy
  std::vector<std::string> _nodal_variables;
  /// Elemental variables to copy
  std::vector<std::string> _elemental_variables;
  /// Block names of the file (ids for unnamed blocks), used to look up values in all of them
  std::vector<SubdomainName> _file_blocks;
  /// Whether the header has been read
  bool _header_read;
};
