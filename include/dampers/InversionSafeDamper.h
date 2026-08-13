#pragma once

#include "GeneralDamper.h"
#include "MooseVariable.h"

#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/tensor_value.h"

class FEProblemBase;
class MooseMesh;
class DisplacedProblem;
class Assembly;

/**
 * Damps the Newton step so that no element inverts on the displaced mesh.
 *
 * MOOSE applies the damped step as x_new = x_old - damping*update, so the deformed configuration at
 * step scale lambda has FE-map Jacobian J(lambda) = J0 - lambda*dJ (note the minus sign), a low-order
 * polynomial in lambda. The base (old) configuration is reconstructed purely from the vectors that
 * the post-check hands us -- x_old = x_ref + (solution + update) -- so the result does not depend on
 * whatever sync state the displaced mesh happens to be in. For each element/quadrature point we find
 * the smallest lambda in (0,1] at which det drops to `min_relative_jacobian` times its current value
 * and return `safety_factor` times the global minimum. This is done purely algebraically from the
 * master-element shape gradients and the nodal update -- it never reinitializes a folded element,
 * so (unlike ElementJacobianDamper) it caps an inverting step instead of failing on it.
 */
class InversionSafeDamper : public GeneralDamper
{
public:
  static InputParameters validParams();

  InversionSafeDamper(const InputParameters & parameters);

  virtual Real computeDamping(const NumericVector<Number> & solution,
                              const NumericVector<Number> & update) override;

protected:
  /// Thread id
  const THREAD_ID _tid;
  /// The FE problem and the (undisplaced/reference) mesh whose map can fold on the displaced config
  FEProblemBase & _fe_problem;
  MooseMesh * _mesh;

  /// The displacement variables
  std::vector<MooseVariable *> _disp_var;
  unsigned int _ndisp;

  /// Allow at most this fraction of the step-to-inversion
  const Real _safety_factor;
  /// Limit before det(J) drops below this fraction of its current value
  const Real _min_relative_jacobian;
  /// Number of lambda samples used to bracket the first crossing in (0,1]
  const unsigned int _n_samples;

  /// libMesh FE used only to obtain the (geometry-independent) master-element shape gradients
  std::unique_ptr<libMesh::FEBase> _fe;
  std::unique_ptr<libMesh::QBase> _qrule;
  const std::vector<std::vector<Real>> * _dphidxi;
  const std::vector<std::vector<Real>> * _dphideta;
  const std::vector<std::vector<Real>> * _dphidzeta;
  unsigned int _dim;
};
