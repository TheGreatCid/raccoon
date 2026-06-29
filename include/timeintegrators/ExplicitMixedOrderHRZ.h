//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#pragma once

#include "ExplicitMixedOrder.h"

/**
 * Variant of MOOSE's ExplicitMixedOrder central-difference integrator that
 * adds the Hinton-Rock-Zienkiewicz (HRZ) mass-matrix lumping scheme.
 *
 * MOOSE's default is row-sum lumping (lumped_i = sum_j M_ij), which produces
 * near-zero or negative diagonal entries for TET10's corner nodes (the
 * consistent mass matrix has negative off-diagonals that overwhelm the
 * small positive diagonal).  The resulting huge 1/m_lumped at corners
 * amplifies high-frequency mass-imbalance modes that central difference
 * cannot damp, and multi-element TET10 meshes go unstable within a few
 * hundred steps.
 *
 * HRZ lumping (Hinton-Rock-Zienkiewicz, 1976) replaces this with:
 *
 *     lumped_i = M_ii * (m_total / sum_j M_jj)
 *     m_total  = sum_ij M_ij      (= total system mass)
 *
 * i.e., take the actual positive diagonal entries and rescale them so their
 * sum equals the total mass.  Guarantees positivity, preserves total mass
 * exactly, and gives a ratio between mid-edge and corner node masses (~4:1
 * for TET10) that's well-conditioned for central difference.
 *
 * Set `mass_lumping_type = hrz` in the [TimeIntegrator] block to enable.
 * Default is `row_sum`, which reproduces the parent class's behavior exactly.
 */
class ExplicitMixedOrderHRZ : public ExplicitMixedOrder
{
public:
  static InputParameters validParams();

  ExplicitMixedOrderHRZ(const InputParameters & parameters);

  virtual void solve() override;

protected:
  /// Lumping method: "row_sum" (parent's behavior) or "hrz".
  const MooseEnum _mass_lumping_type;
};
