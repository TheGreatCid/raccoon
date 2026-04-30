//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"
#include "MaterialProperty.h"
#include "MooseArray.h"
#include "SolutionUserObject.h"
#include "Qp_Mapping.h"
#include "libmesh/dense_matrix.h"
/**
 * This class computes the deformation gradient
 */
class ComputeDeformationGradient : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ComputeDeformationGradient(const InputParameters & parameters);

  void initialSetup() override;

  void computeProperties() override;

  // void computeQpProperties() override;

  void initStatefulProperties(unsigned int n_points) override;

protected:
  // virtual void initQpStatefulProperties() override;

  virtual void displacementIntegrityCheck();

  virtual ADReal computeQpOutOfPlaneGradDisp();

  /// The coordinate system
  const Moose::CoordinateSystemType & _coord_sys;

  /// Coupled displacement variables
  const unsigned int _ndisp;

  /// Displacement variables
  std::vector<const ADVariableValue *> _disp;

  /// Gradient of displacements
  std::vector<const ADVariableGradient *> _grad_disp;

  /// Whether to apply volumetric locking correction
  const bool _volumetric_locking_correction;

  /// The current element volume
  const Real & _current_elem_volume;

  /// The total deformation gradient
  ADMaterialProperty<RankTwoTensor> & _F;
  ADMaterialProperty<RankTwoTensor> & _Fnobar;
  ADMaterialProperty<RankTwoTensor> & _F_store_Fbar;
  const MaterialProperty<RankTwoTensor> & _F_store_Fbar_old;
  // The mechanical deformation gradient (after excluding eigen deformation gradients from the total
  // deformation gradient)
  ADMaterialProperty<RankTwoTensor> & _Fm;

  ADMaterialProperty<RankTwoTensor> & _F_store_noFbar;
  const MaterialProperty<RankTwoTensor> & _F_store_noFbar_old;

  // @{ Eigen deformation gradients
  std::vector<MaterialPropertyName> _Fg_names;
  std::vector<const ADMaterialProperty<RankTwoTensor> *> _Fgs;
  // @}

  ADMaterialProperty<Real> & _weights;

  // is this recovering?
  const bool _recover;

  const SolutionUserObject * _solution_object_ptr;

  std::vector<const ADVariableValue *> _F_recover;

  QpMapping::Element _element = QpMapping::Element::HEX8_3rd;

  unsigned int _qpnum = 0;
  MaterialProperty<Real> & _Frobenius;
  MaterialProperty<Real> & _Jacobian;

  // Polar decomposition output (R and U from F = R*U)
  ADMaterialProperty<RankTwoTensor> & _rotation_tensor;
  const MaterialProperty<RankTwoTensor> & _rotation_tensor_old;
  ADMaterialProperty<RankTwoTensor> & _stretch_tensor;
  /// Volume-corrected stretch tensor: U_fbar = R^T * F_fbar = U * cbrt(<J>_elem / J_qp).
  /// SPD with det = <J>_elem; stored so a downstream restart can recover F_fbar = R * U_fbar
  /// directly (rather than reconstructing F_raw = R * U and re-applying fbar against the
  /// new mesh's _JxW, which integrates over the wrong reference configuration).
  ADMaterialProperty<RankTwoTensor> & _stretch_tensor_fbar;

  // Recovery mode flag
  const bool _recover_from_polar;

  /// When true, store R^(1/2) (rotation by theta/2) instead of R in the output.
  /// Halving the rotation angle keeps it in [0, pi/2), improving conditioning
  /// of the matrix logarithm used by the remapping algorithm.
  const bool _output_half_rotation;

  /// When true, the rotation tensor being read from the recovery file contains R^(1/2)
  /// and must be squared before reconstructing F = R*U.
  /// Set this to match what output_half_rotation_tensor was in the run that produced the file.
  const bool _input_half_rotation;

  /// Use Higham's iterative algorithm for polar decomposition instead of the default
  /// eigendecomposition of C = F^T*F.  Higham operates directly on F (kappa(F) vs kappa(F)^2),
  /// avoids inverting U, requires only 3x3 matrix inverses, and converges quadratically.
  const bool _use_iterative_polar;

  /// Order of the per-element manifold-aware smoothing applied to the recovered R/U/U_fbar
  /// fields at INITIAL.  0 = none (per-QP from solution UO), 1 = linear-in-position fit,
  /// 2 = full quadratic-in-position fit (constant + linear + xx, yy, zz, xy, xz, yz in 3D).
  /// The fit is performed in the log space of each tensor's natural manifold (so(3) for R,
  /// sym(3) for U / U_fbar) and exp'd back.  The motivation is that on cross-mesh restart
  /// the per-QP recovered F field generally does not lie in the range of "I + grad(u)" for
  /// any disp field representable by the new mesh's element type; on TET10_4th especially
  /// this manifests as oscillations on the mid-edge dofs at step 1 as Newton tries to absorb
  /// the per-QP F noise.  Polynomial fitting collapses F_recovered onto the subspace the new
  /// mesh's grad(u) can match.  Quadratic order targets TET10's quadratic disp shape
  /// functions (which yield linear-in-x grad(u) and therefore can match linear-in-x F);
  /// linear-in-x F corresponds to a quadratic-in-x log-tensor field after the polar
  /// decomposition's nonlinearity is unwrapped, hence the quadratic basis on log-space.
  /// The pipeline gracefully drops the order if the element has fewer QPs than the basis.
  unsigned int _recover_smoothing_order;

private:
  const std::unordered_map<int, int> * _lookup;

  /// Compute the principal square root of a rotation matrix R in SO(3).
  /// R_half_old is the stored rotation_tensor from the previous timestep; its skew
  /// part is used to detect and correct axis-sign flips that occur when the physical
  /// rotation crosses π, maintaining component continuity across timesteps.
  static RankTwoTensor computeHalfRotation(const RankTwoTensor & R,
                                            const RankTwoTensor & R_half_old);

  /// Higham iterative polar decomposition: F = R * U where R in SO(3) and U is symmetric PD.
  /// X_{k+1} = (X_k + X_k^{-T}) / 2, converges quadratically to R.
  /// More numerically stable than getRUDecompositionRotation for ill-conditioned F.
  void polarDecompositionIterative(const RankTwoTensor & F, RankTwoTensor & R, RankTwoTensor & U);

  /// Matrix log of an SPD 3x3 tensor M, via symmetric eigendecomposition:
  /// log(M) = sum_k log(lambda_k) v_k v_k^T.  Errors out if any eigenvalue is non-positive.
  static RankTwoTensor logSPD(const RankTwoTensor & M);

  /// Matrix exp of a symmetric 3x3 tensor S, via symmetric eigendecomposition:
  /// exp(S) = sum_k exp(lambda_k) v_k v_k^T.
  static RankTwoTensor expSym(const RankTwoTensor & S);

  /// Matrix log of a rotation R in SO(3); returns the skew-symmetric log in so(3).
  /// Uses Rodrigues' inverse for general theta and the symmetric-part axis trick near theta = pi.
  static RankTwoTensor logSO3(const RankTwoTensor & R);

  /// Matrix exp of a skew-symmetric 3x3 tensor K; returns the rotation in SO(3).
  /// Uses Rodrigues with theta = ||K|| / sqrt(2); falls back to Taylor for small theta.
  static RankTwoTensor expSO3(const RankTwoTensor & K);

  /// Fit each (i, j) component of `log_tensors[qp]` against a linear-in-position
  /// polynomial (constant + ndisp linear terms) over the element's QPs and overwrite
  /// the inputs with the fit values.  `qp_basis[qp]` is the precomputed basis vector
  /// at QP `qp` of size `nb`.  `lu` is the precomputed LU of the normal-equations
  /// matrix sum_qp basis_qp basis_qp^T.  If `enforce` is "skew" or "sym", each
  /// fit tensor is projected onto the corresponding subspace before returning.
  enum class LogProjection { None, Skew, Sym };
  static void fitLogField(std::vector<RankTwoTensor> & log_tensors,
                          const std::vector<std::vector<Real>> & qp_basis,
                          const DenseMatrix<Real> & normal_lu,
                          unsigned int nb,
                          LogProjection projection);
};
