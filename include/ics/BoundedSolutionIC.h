#pragma once

#include "SolutionIC.h"

/**
 * A SolutionIC that clamps the recovered value into [lower_bound, upper_bound].
 *
 * On restart the shape-evaluation remapper can return a phase-field value slightly
 * outside its physical range (e.g. a small negative d, or an overshoot above 1). Seeding
 * the initial condition -- and, for a VI-bounded variable, the recovered lower bound --
 * with such a value is undesirable. This reuses SolutionIC's recovery logic verbatim and
 * only clamps the returned scalar.
 */
class BoundedSolutionIC : public SolutionIC
{
public:
  static InputParameters validParams();

  BoundedSolutionIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

protected:
  /// Recovered values below this are clamped up to it
  const Real _lower_bound;
  /// Recovered values above this are clamped down to it
  const Real _upper_bound;
};
