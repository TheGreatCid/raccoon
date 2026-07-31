#pragma once

#include "SolutionUserObject.h"
#include "Qp_Mapping.h"

class SolutionUserObjectQP : public SolutionUserObject
{
public:
  static InputParameters validParams();
  SolutionUserObjectQP(const InputParameters & parameters);

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
  const std::unordered_map<int, int> * _lookup;
};
