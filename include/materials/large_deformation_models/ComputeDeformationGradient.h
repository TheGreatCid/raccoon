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

  /// Impose the out-of-plane component of the displacement-gradient tensor A (before the identity is
  /// added, so F_zz = 1 + A(2,2)). The base applies it only for RZ (axisymmetric hoop stretch);
  /// derived plane-stress models override this to set A(2,2) from a coupled out-of-plane strain.
  virtual void applyOutOfPlaneGradDisp(ADRankTwoTensor & A);

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
  /// det(F_raw) per QP -- declared AD so the RecoverVariables/exodus dump
  /// pipeline (which wires ADMaterialRealAux) can pick it up, and so the
  /// recovery side reading Jacobian back via SolutionReal has the AD
  /// version it expects.
  ADMaterialProperty<Real> & _Jacobian;

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

  /// When recovering with volumetric-locking correction active, controls how F_bar is
  /// constructed on the restart side from the recovered fields:
  ///   false (default): F_bar = R * U_fbar, where U_fbar is the per-QP F-bar-corrected
  ///                    stretch tensor read directly from the recovery file (approach A).
  ///                    This is what the rest of the current pipeline does and is the
  ///                    behaviour you get if the recovery file already contains
  ///                    `stretch_tensor_fbar`.
  ///   true           : F_bar is rebuilt on the new mesh by applying the original
  ///                    cfb539fe1-era F-bar averaging operator to the recovered raw U:
  ///                       F_raw[qp] = R[qp] * U[qp]
  ///                       J_avg     = (1/V) * sum_qp det(F_raw[qp]) * JxW * coord
  ///                       F_bar[qp] = F_raw[qp] * cbrt(J_avg / det(F_raw[qp]))
  ///                    (approach B).  Useful when the recovery file does NOT carry
  ///                    `stretch_tensor_fbar` and the F-bar correction has to be
  ///                    re-derived from the raw stretch tensor.
  /// Active only when both `recover = true` and `volumetric_locking_correction = true`.
  const bool _recover_apply_fbar_to_U;

  /// Approach C: apply F-bar once to the EXACT total F (= F_inc * F_dump_raw),
  /// using a change-of-variables-corrected J_avg over the original undeformed
  /// volume.  Required for explicit dynamics (no Newton to absorb the
  /// multiplicative F-bar composition error of approach A/B).  Implies raw U
  /// recovery from the dump; mutually exclusive with `_recover_apply_fbar_to_U`.
  const bool _recover_apply_fbar_to_total;

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
};
