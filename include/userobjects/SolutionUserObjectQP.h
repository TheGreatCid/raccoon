//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "SolutionUserObject.h"
#include "Qp_Mapping.h"
#include "RankTwoTensor.h"

/**
 * A SolutionUserObject for recovery files, in which material properties are stored as one
 * elemental variable per quadrature point (see RecoverVariablesAction, which writes them). The
 * variables of the listed properties are read, and can be evaluated at the quadrature points of an
 * element of this simulation. The element is matched to the file by its centroid.
 */
class SolutionUserObjectQP : public SolutionUserObject
{
public:
  static InputParameters validParams();
  SolutionUserObjectQP(const InputParameters & parameters);

  /// A scalar property of the file at a quadrature point (0-based, MOOSE numbering) of an element
  Real qpValue(Real t, const Elem & elem, const std::string & property, unsigned int qp) const;

  /// A tensor property of the file at a quadrature point (0-based, MOOSE numbering) of an element
  RankTwoTensor
  qpTensor(Real t, const Elem & elem, const std::string & property, unsigned int qp) const;

  /// The variable holding a scalar property / a tensor component at a quadrature point (0-based,
  /// MOOSE numbering)
  std::string qpVariableName(const std::string & property, unsigned int qp) const;
  std::string qpVariableName(const std::string & property,
                             unsigned int i,
                             unsigned int j,
                             unsigned int qp) const;

protected:
  /// Releases the large resident payload after the one-time initial recovery reads
  virtual void timestepSetup() override;

  std::vector<MaterialName> _tensor_materials;
  std::vector<MaterialName> _materials;

  QpMapping::Element _element;

  unsigned int _qpnum;

  /// If true, release the loaded solution/mesh payload after the one-time initial
  /// recovery reads. Only safe for the recovery pattern (this object is queried
  /// solely during stateful-property initialization); incompatible with time
  /// interpolation. Do not enable when the object is queried after INITIAL.
  const bool _free_after_initial;

  /// Guard so the payload is released only once
  bool _payload_freed;

private:
  /// Map from the MOOSE quadrature point numbering to the file's
  const std::unordered_map<int, int> * _lookup;
};

