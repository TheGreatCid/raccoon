//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "Material.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"

class SolutionUserObjectQP;

/**
 * This class computes the deformation gradient, with an optional F-bar volumetric locking
 * correction and removal of eigen deformation gradients. It also outputs the polar decomposition
 * F = R * U of the uncorrected deformation gradient, which is what a later run recovers from.
 *
 * With recover = true the mesh is the deformed configuration of an earlier run, the displacements
 * are incremental, and the total deformation gradient is F_inc * F_rec, where F_rec is read at
 * each quadrature point from a SolutionUserObjectQP (as R and U, or as F). When the F-bar
 * correction is active there are three ways to recover the corrected deformation gradient:
 *
 *   - default: read the corrected stretch (stretch_tensor_fbar) and use F_bar_rec = R * U_bar.
 *     The total is F_bar(F_inc) * F_bar_rec.
 *   - recover_apply_fbar_to_U: read the raw stretch (stretch_tensor) and rebuild F_bar_rec by
 *     averaging det(R * U) over the element. The total is F_bar(F_inc) * F_bar_rec.
 *   - recover_apply_fbar_to_total: read the raw stretch and apply F-bar once to the total raw
 *     deformation gradient F_inc * F_rec, averaging over the original undeformed volume.
 */
class ComputeDeformationGradient : public Material, public BaseNameInterface
{
public:
  static InputParameters validParams();

  ComputeDeformationGradient(const InputParameters & parameters);

  void initialSetup() override;

  void computeProperties() override;

  void initStatefulProperties(unsigned int n_points) override;

protected:
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
  /// Not read. Requesting the old value makes the recovered F stateful, so the value seeded in
  /// initStatefulProperties() is kept across time steps.
  const MaterialProperty<RankTwoTensor> & _F_store_Fbar_old;
  // The mechanical deformation gradient (after excluding eigen deformation gradients from the total
  // deformation gradient)
  ADMaterialProperty<RankTwoTensor> & _Fm;

  ADMaterialProperty<RankTwoTensor> & _F_store_noFbar;
  /// Not read; makes the recovered F stateful (see _F_store_Fbar_old)
  const MaterialProperty<RankTwoTensor> & _F_store_noFbar_old;

  // @{ Eigen deformation gradients
  std::vector<MaterialPropertyName> _Fg_names;
  std::vector<const ADMaterialProperty<RankTwoTensor> *> _Fgs;
  // @}

  /// Whether the deformation of an earlier run is recovered
  const bool _recover;

  /// The solution the recovered deformation is read from
  const SolutionUserObjectQP * _solution_object_ptr;

  /// Components of an externally provided deformation gradient to recover with
  std::vector<const ADVariableValue *> _F_recover;

  /// Frobenius norm of F - I
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

  /// Rebuild the recovered F-bar from the raw recovered stretch (see the class description)
  const bool _recover_apply_fbar_to_U;

  /// Apply F-bar once to the total raw deformation gradient (see the class description)
  const bool _recover_apply_fbar_to_total;

  /// The volumetric_locking_correction parameter as given. _volumetric_locking_correction is
  /// false for boundary materials, but which stretch tensor the recovery file holds does not
  /// depend on where this material is evaluated.
  const bool _recover_with_fbar;

private:
  /**
   * Read the recovered rotation and stretches at every quadrature point of the current element.
   * @param R       the rotation
   * @param R_file  the rotation as stored in the file (its square root with input_half_rotation)
   * @param U       the raw stretch (the corrected stretch if the file only holds that)
   * @param U_fbar  the F-bar corrected stretch (the raw stretch if the file only holds that)
   */
  void readRecoveredState(std::vector<RankTwoTensor> & R,
                          std::vector<RankTwoTensor> & R_file,
                          std::vector<RankTwoTensor> & U,
                          std::vector<RankTwoTensor> & U_fbar) const;

  /// Whether the recovery needs the raw / the F-bar corrected stretch from the file
  bool needRawStretch() const;
  bool needFbarStretch() const;

  /// Warn if R is not a rotation or R * U does not reproduce F
  void checkPolarDecomposition(const RankTwoTensor & R,
                               const RankTwoTensor & U,
                               const RankTwoTensor & F) const;

  /// Compute the principal square root of a rotation matrix R in SO(3).
  /// R_half_old is the stored rotation_tensor from the previous timestep; its skew
  /// part is used to detect and correct axis-sign flips that occur when the physical
  /// rotation crosses π, maintaining component continuity across timesteps.
  static RankTwoTensor computeHalfRotation(const RankTwoTensor & R,
                                            const RankTwoTensor & R_half_old);

  /// Higham iterative polar decomposition: F = R * U where R in SO(3) and U is symmetric PD.
  /// X_{k+1} = (X_k + X_k^{-T}) / 2, converges quadratically to R.
  /// More numerically stable than getRUDecompositionRotation for ill-conditioned F.
  void
  polarDecompositionIterative(const RankTwoTensor & F, RankTwoTensor & R, RankTwoTensor & U) const;
};
