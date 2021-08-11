//* David's development

#pragma once

#include "ADKernelValue.h"
#include "BaseNameInterface.h"
#include "ADRankTwoTensorForward.h"

class ADgdeg_T : public ADKernelValue
{
public:
  static InputParameters validParams();
  ADgdeg_T(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;
  virtual ADReal rexp_calc();
  virtual ADReal rtanh_calc();

  const Real _Q;
  const ADMaterialProperty<Real> & _tau_bar;
  const ADMaterialProperty<Real> & _ep;
  const Real _del;
  const Real _kappa;
  const Real _T0;
  const Real _ep0;
  const Real _ep0_dot;
  const Real _sigma0;
  const Real _sigma_eff;
  const Real _k;
};
