//* This file is part of the RACCOON application
//* being developed at Dolbow lab at Duke University
//* http://dolbow.pratt.duke.edu

#include "ADRankTwoTensorForward.h"
#include "ADReal.h"
#include "EigenADReal.h"
#include "LargeDeformationJ2Plasticity.h"
#include "MooseError.h"
#include "RaccoonUtils.h"
#include <cmath>
#include <string>

registerMooseObject("raccoonApp", LargeDeformationJ2Plasticity);

InputParameters
LargeDeformationJ2Plasticity::validParams()
{
  InputParameters params = LargeDeformationPlasticityModel::validParams();
  params.addClassDescription("Large deformation $J_2$ plasticity. The exponential constitutive "
                             "update is used to update the plastic deformation.");
  params.addParam<bool>("recover", false, "do you want to recover");
  params.addParam<bool>("coupled_temp_solve",
                        false,
                        "Solve for plastic strain increment and adiabatic temperature rise "
                        "simultaneously using a local 2x2 Newton, rather than using the lagged "
                        "temperature from the hardening model.");
  params.addParam<MaterialPropertyName>(
      "density", "rho", "Mass density material property (required when coupled_temp_solve=true).");
  params.addParam<Real>(
      "specific_heat",
      0.0,
      "Specific heat capacity c_v [J/(kg*K)] (required when coupled_temp_solve=true).");
  return params;
}

LargeDeformationJ2Plasticity::LargeDeformationJ2Plasticity(const InputParameters & parameters)
  : LargeDeformationPlasticityModel(parameters),
    _phi(declareADProperty<Real>("phi")),
    _flowstress(declareADProperty<Real>("flowstress")),
    _visflowstress(declareADProperty<Real>("visflowstress")),
    _coupled_temp_solve(getParam<bool>("coupled_temp_solve")),
    _rho(nullptr),
    _cv(getParam<Real>("specific_heat")),
    _dep_2x2_guess(0.0),
    _use_2x2_guess(false)
{
  _check_range = true;
  if (_coupled_temp_solve)
  {
    if (_cv == 0.0)
      mooseWarning(name(),
                   ": coupled_temp_solve=true but specific_heat=0. "
                   "Temperature rise will be zero. Set specific_heat to enable thermal coupling.");
    _rho = &getADMaterialProperty<Real>(getParam<MaterialPropertyName>("density"));
  }
}

ADReal
LargeDeformationJ2Plasticity::initialGuess(const ADReal & effective_trial_stress)
{
  if (_use_2x2_guess)
    return ADReal(_dep_2x2_guess);
  return LargeDeformationPlasticityModel::initialGuess(effective_trial_stress);
}

void
LargeDeformationJ2Plasticity::updateState(ADRankTwoTensor & stress, ADRankTwoTensor & Fe)
{

  ADRankTwoTensor curr_Fp;

  // populate curr_FP and _ep_old_store
  //_ep_old store is a material property so that it can be used in the residual calculation
  using std::sqrt;

  // First assume no plastic increment
  ADReal delta_ep = 0;
  Fe = Fe * _Fp_old[_qp].inverse();

  stress = _elasticity_model->computeMandelStress(Fe);
  // Compute the flow direction following the Prandtl-Reuss flow rule.
  // We guard against zero denominator.
  ADRankTwoTensor stress_dev = stress.deviatoric();
  ADReal stress_dev_norm = stress_dev.doubleContraction(stress_dev);
  if (MooseUtils::absoluteFuzzyEqual(stress_dev_norm, 0))
    stress_dev_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
  stress_dev_norm = sqrt(1.5 * stress_dev_norm);
  _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;

  _phi[_qp] = computeResidual(stress_dev_norm, delta_ep);

  if (_coupled_temp_solve)
  {
    if (_phi[_qp] > 0)
    {
      // ---- Coupled 2x2 local Newton for (delta_ep, delta_T) ----
      // R1 = yield residual = 0
      // R2 = rho*cv*delta_T - (W_d + TC) * delta_ep = 0  (adiabatic heat balance)
      const Real T_old = _hardening_model->getQpTemperatureOld();
      const Real rho_cv = MetaPhysicL::raw_value((*_rho)[_qp]) * _cv;
      const Real eff_trial = MetaPhysicL::raw_value(stress_dev_norm);
      const Real R1_ref = std::max(eff_trial, 1.0);

      // Evaluate initial guess using the old timestep temperature (safer starting point)
      _hardening_model->setLocalTemperature(T_old);
      Real dep_raw = MetaPhysicL::raw_value(_hardening_model->initialGuess(stress_dev_norm));
      _hardening_model->clearLocalTemperature();

      // Clamp the scalar initial guess to a sensible value (debug / safety; tune as needed)
      const Real dep_max = 1.0; // absolute max plastic strain increment allowed as an initial guess
      const Real dep_safe_init = std::max(std::min(dep_raw, dep_max), 1e-12);
      Real dep = dep_safe_init;
      Real dT = 0.0;

      // Evaluates both residuals at a given (dep, dT), injecting the temperature into the
      // hardening model before each evaluation. Also returns Y, W_d, TC for use in the
      // analytic Jacobian.
      auto evalR = [&](Real dep_in,
                       Real dT_in,
                       Real & R1,
                       Real & R2,
                       Real & Y_out,
                       Real & W_d_out,
                       Real & TC_out)
      {
        Real T = T_old + dT_in;
        Real ep = MetaPhysicL::raw_value(_ep_old[_qp]) + dep_in;
        if (ep <= 0.0)
          ep = 1e-20;
        const Real dep_safe = (dep_in <= 0.0) ? 1e-20 : dep_in;

        _hardening_model->setLocalTemperature(T);

        const Real elas_corr = MetaPhysicL::raw_value(
            _elasticity_model->computeMandelStress(dep_safe * _Np[_qp], /*plasticity_update=*/true)
                .doubleContraction(_Np[_qp]));
        const Real Y = MetaPhysicL::raw_value(_hardening_model->plasticEnergy(ADReal(ep), 1));
        const Real W_d = MetaPhysicL::raw_value(
            _hardening_model->plasticDissipation(ADReal(dep_safe), ADReal(ep), 1));
        const Real TC = MetaPhysicL::raw_value(_hardening_model->thermalConjugate(ADReal(ep)));

        R1 = eff_trial - elas_corr - Y - W_d;
        R2 = rho_cv * dT_in - (W_d + TC) * dep_safe;
        Y_out = Y;
        W_d_out = W_d;
        TC_out = TC;
      };

      const unsigned int max_iter = 200;
      const Real tol = 1e-12;

      // Jacobian at the converged point — populated inside the loop, used after for the
      // consistent tangent correction applied to the AD return-mapping derivatives.
      Real J00_f = 0, J01_f = 0, J10_f = 0, J11_f = 0, det_f = 1;

      for (unsigned int iter = 0; iter < max_iter; ++iter)
      {
        Real R1, R2, Y, W_d, TC;
        evalR(dep, dT, R1, R2, Y, W_d, TC);

        // Relative convergence: normalize R1 by stress scale, R2 by rho*cv (1 K rise)
        const Real R2_ref = std::max(rho_cv, 1.0);
        const Real rel_norm =
            std::sqrt((R1 / R1_ref) * (R1 / R1_ref) + (R2 / R2_ref) * (R2 / R2_ref));

        // Analytic Jacobian — evaluated unconditionally so the converged J is available
        // outside the loop for the global tangent correction.
        const Real dep_safe_j = (dep <= 0.0) ? 1e-20 : dep;
        const Real ep_j = std::max(MetaPhysicL::raw_value(_ep_old[_qp]) + dep_safe_j, 1e-20);
        const Real T_curr = T_old + dT;
        // xi = d(TD)/dT / TD; used for dY/dT = Y*xi and dW_d/dT = W_d*xi
        const Real xi = _hardening_model->temperatureDependenceLogDerivative(T_curr);
        // J[0][0] = dR1/d(dep): reuse the existing computeDerivative (analytic, exact)
        J00_f = MetaPhysicL::raw_value(computeDerivative(stress_dev_norm, ADReal(dep_safe_j)));
        // J[0][1] = dR1/d(dT) = -(dY/dT + dW_d/dT) = -(Y + W_d)*xi
        J01_f = -(Y + W_d) * xi;
        // J[1][0] = dR2/d(dep) = -(W_d + TC) - dep * d(W_d+TC)/d(dep)
        // Use the full total derivative d(W_d)/d(dep) = plasticDissipation(dep,ep,2),
        // which includes both the direct dep term (rate) and the ep = ep_old+dep term
        // (hardening). Also include dep * dTC/dep via the ep channel.
        const Real dWd_dep_total =
            dep_safe_j * MetaPhysicL::raw_value(_hardening_model->plasticDissipation(
                             ADReal(dep_safe_j), ADReal(ep_j), 2));
        const Real dTC_dep = _hardening_model->thermalConjugateEpDerivative(ep_j);
        J10_f = -(W_d + TC + dWd_dep_total + dep_safe_j * dTC_dep);
        // J[1][1] = dR2/d(dT) = rho_cv - dep*(dW_d/dT + dTC/dT)
        const Real dTC_dT = _hardening_model->thermalConjugateTemperatureDerivative(ep_j);
        J11_f = rho_cv - dep_safe_j * (W_d * xi + dTC_dT);
        det_f = J00_f * J11_f - J01_f * J10_f;

        if (rel_norm < tol)
          break;

        if (iter == max_iter - 1)
          mooseError(name(),
                     ": coupled thermo-plastic local Newton did not converge at qp ",
                     _qp,
                     ". Final relative residual norm: ",
                     rel_norm,
                     ". R1=",
                     R1,
                     " R2=",
                     R2);

        if (std::abs(det_f) < 1e-30 * std::max(std::abs(J00_f * J11_f), 1.0))
          mooseError(
              name(), ": singular local Jacobian in coupled thermo-plastic solve at qp ", _qp);

        // Newton update: dx = -J^{-1} * R
        // Limit the dep step to at most 50% relative change per iteration.
        // This prevents the 1/dep singularity in the rate-dependent log term from
        // making J[0][0] astronomically large when dep→0, which would reduce the
        // Newton step to near zero and stall convergence.
        const Real step_dep_full = -(J11_f * R1 - J01_f * R2) / det_f;
        const Real step_dT_full = -(J00_f * R2 - J10_f * R1) / det_f;
        const Real dep_step_limit = 0.5 * std::max(dep, 1e-14);
        const Real step_dep = std::max(std::min(step_dep_full, dep_step_limit), -dep_step_limit);
        dep = std::max(dep + step_dep, 1e-14);
        dT += step_dT_full;
      }

      // Fix temperature at converged value so stress and _heat are evaluated at T_old + dT.
      // With the new getQpT() offset trick, this also threads ∂/∂T_global through all
      // subsequent material evaluations (plasticDissipation, thermalConjugate, etc.).
      _hardening_model->setLocalTemperature(T_old + dT);

      // Construct delta_ep as an ADReal with the exact 2×2 consistent tangent.
      //
      // From the implicit function theorem on the 2×2 system:
      //   ∂dep/∂σ_eff  = -J11 / det   (adiabatic correction to the mechanical tangent)
      //   ∂dep/∂T_g    =  J01 * ρcv / det  (sensitivity to the global temperature DOF)
      //
      // Expressed as a first-order expansion around the converged (dep, T_g) values,
      // using zero-value/nonzero-derivative ADReal adjustments:
      //   sig_adj carries ∂σ_eff/∂u_DOFs with value 0
      //   T_adj   carries ∂T_g/∂T_DOF  = 1 with value 0
      //
      // This replaces the previous returnMappingSolve + scalar cf correction, which
      // applied the same factor to ALL AD channels and gave the wrong sign for ∂dep/∂T_g.
      const ADReal T_ad = _hardening_model->getQpTAD();
      const ADReal sig_adj = stress_dev_norm - ADReal(MetaPhysicL::raw_value(stress_dev_norm));
      const ADReal T_adj = T_ad - ADReal(MetaPhysicL::raw_value(T_ad));

      if (std::abs(det_f) > 1e-30 * std::max(std::abs(J00_f * J11_f), 1.0))
      {
        delta_ep = ADReal(dep)
                   + (-J11_f / det_f) * sig_adj
                   + (J01_f * rho_cv / det_f) * T_adj;
      }
      else
      {
        // Singular or near-singular 2×2 Jacobian: fall back to uncoupled 1D tangent
        // (∂dep/∂σ_eff = -1/J00, no T-channel correction).
        delta_ep = ADReal(dep) + (-1.0 / J00_f) * sig_adj;
      }
    }
  }
  else
  {
    // ensure hardening uses lagged temperature explicitly
    const Real T_old = _hardening_model->getQpTemperatureOld();
    _hardening_model->setLocalTemperature(T_old);

    if (_phi[_qp] > 0)
      returnMappingSolve(stress_dev_norm, delta_ep, _console);
  }

  _ep[_qp] = _ep_old[_qp] + delta_ep;

  // Avoiding /0 issues for rate dependent models
  _ep[_qp] = (_ep[_qp] == 0) ? 1e-20 : _ep[_qp];

  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  Fe = Fe * delta_Fp.inverse();
  stress = _elasticity_model->computeCauchyStress(Fe);

  _hardening_model->plasticEnergy(_ep[_qp]);
  _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 0);

  // Avoiding NaN issues for rate depedent models
  if (_t_step > 0)
  {
    _heat[_qp] = _hardening_model->plasticDissipation(delta_ep, _ep[_qp], 1) * delta_ep / _dt;

    _heat[_qp] += _hardening_model->thermalConjugate(_ep[_qp]) * delta_ep / _dt;
  }
  else
    _heat[_qp] = 0;

  // Release any local temperature override so the next call starts clean.
  // This is done here (not inside the branches) so that plasticEnergy,
  // plasticDissipation, and _heat above are all evaluated at the same
  // temperature that was used in the return mapping:
  //   coupled path:     T_old + dT  (set by setLocalTemperature before returnMappingSolve)
  //   non-coupled path: T_old       (set by setLocalTemperature at the top of the else branch)
  _hardening_model->clearLocalTemperature();
}

Real
LargeDeformationJ2Plasticity::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                       const ADReal & delta_ep)
{
  return raw_value(effective_trial_stress - _elasticity_model
                                                ->computeMandelStress(delta_ep * _Np[_qp],
                                                                      /*plasticity_update = */ true)
                                                .doubleContraction(_Np[_qp]));
}

ADReal
LargeDeformationJ2Plasticity::computeResidual(const ADReal & effective_trial_stress,
                                              const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;

  // Avoiding /0 errors for rate depedent models
  ep = (ep == 0) ? 1e-20 : ep;

  return effective_trial_stress -
         _elasticity_model
             ->computeMandelStress(delta_ep * _Np[_qp],
                                   /*plasticity_update = */ true)
             .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(ep, 1) -
         _hardening_model->plasticDissipation(delta_ep, ep, 1);
}

ADReal
LargeDeformationJ2Plasticity::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                const ADReal & delta_ep)
{
  ADReal ep = _ep_old[_qp] + delta_ep;
  if (ep == 0)
  {
    ep = 1e-20;
  }
  return -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
              .doubleContraction(_Np[_qp]) -
         _hardening_model->plasticEnergy(ep, 2) -
         _hardening_model->plasticDissipation(delta_ep, ep, 2);
}
