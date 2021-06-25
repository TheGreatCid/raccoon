#include "ComputeLargeDeformationXMATStress.h"

registerMooseObject("raccoonApp", ComputeLargeDeformationXMATStress);

InputParameters
ComputeLargeDeformationXMATStress::validParams()
{
  InputParameters params = Material::validParams();
  params += BaseNameInterface::validParams();
  params += ADSingleVariableReturnMappingSolution::validParams();

  params.addClassDescription("Stress calculator given an elasticity model, a plasticity model and "
                             "a viscoelasticity model. Large deformation is assumed.");

  // Elasticity
  params.addRequiredParam<MaterialPropertyName>("bulk_modulus", "The bulk modulus $K$");
  params.addRequiredParam<MaterialPropertyName>("shear_modulus", "The shear modulus $G$");

  // Fracture
  params.addRequiredCoupledVar("phase_field", "Name of the phase-field (damage) variable");
  params.addParam<MaterialPropertyName>(
      "strain_energy_density",
      "psie",
      "Name of the strain energy density computed by this material model");
  params.addParam<MaterialPropertyName>("degradation_function", "g", "The degradation function");

  // Creep
  params.addRequiredParam<Real>("coefficient",
                                "The coefficient $A$ in the update formula for creep rate");
  params.addRequiredParam<Real>("exponent",
                                "The exponent $n$ in the update formula for creep rate");
  params.addRequiredParam<MaterialPropertyName>("reference_stress",
                                                "The reference stress $\\sigma_0$");
  params.addRequiredParam<MaterialPropertyName>("arrhenius_coefficient",
                                                "The arrhenius coefficient");
  params.addRequiredParam<Real>(
      "eps", "A small number to help this perfect plasticity model to converge.");
  params.addParam<MaterialPropertyName>(
      "plastic_energy_density",
      "psip",
      "Name of the plastic energy density computed by this material model");
  params.addParam<MaterialPropertyName>(
      "plastic_degradation_function", "gp", "The degradation function");

  params.suppressParameter<bool>("use_displaced_mesh");

  // David's Stuff-------------------
  params.addParam<MaterialPropertyName>("internal_dissipation_energy",
                                        "_del_delt",
                                        "The internal Dissipation energy for heat conduction");
  params.addRequiredParam<MaterialName>("elasticity_model",
                                        "Name of the elastic stress-strain constitutive model");
  //-- ---------------------------------
  return params;
}

ComputeLargeDeformationXMATStress::ComputeLargeDeformationXMATStress(
    const InputParameters & parameters)
  : Material(parameters),
    BaseNameInterface(parameters),
    ADSingleVariableReturnMappingSolution(parameters),
    DerivativeMaterialPropertyNameInterface(),
    // David's Stuff
    _del_delt(declareADProperty<Real>("internal_dissipation_energy")),

    _stress(declareADProperty<RankTwoTensor>(prependBaseName("stress"))),
    //--------
    _Fm(getADMaterialProperty<RankTwoTensor>(prependBaseName("mechanical_deformation_gradient"))),

    ///////////////////////////////
    // Elasticity
    ///////////////////////////////
    _Fe(declareADProperty<RankTwoTensor>(prependBaseName("elastic_deformation_gradient"))),
    _K(getADMaterialPropertyByName<Real>(prependBaseName("bulk_modulus", true))),
    _G(getADMaterialPropertyByName<Real>(prependBaseName("shear_modulus", true))),

    ///////////////////////////////
    // Fracture
    ///////////////////////////////
    _d_name(getVar("phase_field", 0)->name()),
    // The strain energy density and its derivatives
    _psie_name(prependBaseName("strain_energy_density", true)),
    _psie(declareADProperty<Real>(_psie_name)),
    _psie_active(declareADProperty<Real>(_psie_name + "_active")),
    _dpsie_dd(declareADProperty<Real>(derivativePropertyName(_psie_name, {_d_name}))),
    // The degradation function and its derivatives
    _g_name(prependBaseName("degradation_function", true)),
    _g(getADMaterialProperty<Real>(_g_name)),
    _dg_dd(getADMaterialProperty<Real>(derivativePropertyName(_g_name, {_d_name}))),

    ///////////////////////////////
    // Plasticity
    ///////////////////////////////
    _Fp(declareADProperty<RankTwoTensor>(prependBaseName("plastic_deformation_gradient"))),
    _Fp_old(getMaterialPropertyOldByName<RankTwoTensor>(
        prependBaseName("plastic_deformation_gradient"))),
    _ep(declareADProperty<Real>(prependBaseName("effective_plastic_strain"))),
    _ep_old(getMaterialPropertyOldByName<Real>(prependBaseName("effective_plastic_strain"))),
    _Np(declareADProperty<RankTwoTensor>(prependBaseName("flow_direction"))),

    ///////////////////////////////
    // Creep
    ///////////////////////////////
    _coefficient(getParam<Real>("coefficient")),
    _exponent(getParam<Real>("exponent")),
    _sigma_0(getADMaterialProperty<Real>(prependBaseName("reference_stress", true))),
    _arrhenius_coef(getADMaterialProperty<Real>(prependBaseName("arrhenius_coefficient", true))),
    _eps(getParam<Real>("eps")),
    // The strain energy density and its derivatives
    _psip_name(prependBaseName("plastic_energy_density", true)),
    _psip(declareADProperty<Real>(_psip_name)),
    _psip_active(declareADProperty<Real>(_psip_name + "_active")),
    _dpsip_dd(declareADProperty<Real>(derivativePropertyName(_psip_name, {_d_name}))),
    // The degradation function and its derivatives
    _gp_name(prependBaseName("plastic_degradation_function", true)),
    _gp(getADMaterialProperty<Real>(_gp_name)),
    _dgp_dd(getADMaterialProperty<Real>(derivativePropertyName(_gp_name, {_d_name})))

{
  if (getParam<bool>("use_displaced_mesh"))
    mooseError("The stress calculator needs to run on the undisplaced mesh.");
}

void
ComputeLargeDeformationXMATStress::initQpStatefulProperties()
{
  _stress[_qp].zero();
  _Fp[_qp].setToIdentity();
  _ep[_qp] = 0;
}

void
ComputeLargeDeformationXMATStress::computeQpProperties()
{
  initialSetup();
  _Fe[_qp] = _Fm[_qp];
  updateState();
}

void
ComputeLargeDeformationXMATStress::updateState()
{
  // First assume no plastic increment
  ADReal delta_ep = 0;
  _Fe[_qp] = _Fe[_qp] * _Fp_old[_qp].inverse();
  //_stress[_qp] = MandelStress(_Fe[_qp]);
  _stress[_qp] = _elasticity_model->computeMandelStress(_Fe[_qp]);
  // Compute the flow direction following the Prandtl-Reuss flow rule.
  // We guard against zero denominator.
  ADRankTwoTensor stress_dev = _stress[_qp].deviatoric();
  ADReal stress_dev_norm = stress_dev.doubleContraction(stress_dev);
  if (MooseUtils::absoluteFuzzyEqual(stress_dev_norm, 0))
    stress_dev_norm.value() = libMesh::TOLERANCE * libMesh::TOLERANCE;
  stress_dev_norm = std::sqrt(1.5 * stress_dev_norm);
  _Np[_qp] = 1.5 * stress_dev / stress_dev_norm;

  // Return mapping
  ADReal phi = computeResidual(stress_dev_norm, delta_ep);
  if (phi > 0)
    returnMappingSolve(stress_dev_norm, delta_ep, _console);
  _ep[_qp] = _ep_old[_qp] + delta_ep;
  ADRankTwoTensor delta_Fp = RaccoonUtils::exp(delta_ep * _Np[_qp]);
  _Fp[_qp] = delta_Fp * _Fp_old[_qp];

  // Update stress
  _Fe[_qp] = _Fe[_qp] * delta_Fp.inverse();
  _stress[_qp] = CauchyStress(_Fe[_qp]);
  plasticEnergy(_ep[_qp]);
  computeDelValues(delta_ep);
}

ADRankTwoTensor
ComputeLargeDeformationXMATStress::CauchyStress(const ADRankTwoTensor & Fe)
{
  return Fe.inverse().transpose() * _elasticity_model->computeMandelStress(_Fe[_qp]) *
         Fe.transpose() / Fe.det();
}

Real
ComputeLargeDeformationXMATStress::computeReferenceResidual(const ADReal & effective_trial_stress,
                                                            const ADReal & delta_ep)
{
  return raw_value(
      effective_trial_stress -
      _elasticity_model->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
          .doubleContraction(_Np[_qp]));
}

ADReal
ComputeLargeDeformationXMATStress::computeResidual(const ADReal & effective_trial_stress,
                                                   const ADReal & delta_ep)
{
  const ADReal stress_delta =
      effective_trial_stress -
      _elasticity_model->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
          .doubleContraction(_Np[_qp]);
  const ADReal yield_stress = plasticEnergy(_ep_old[_qp] + delta_ep, 1);
  const ADReal creep_rate = _coefficient * std::pow(stress_delta / yield_stress, _exponent);
  return creep_rate * _dt - delta_ep;
}

ADReal
ComputeLargeDeformationXMATStress::computeDerivative(const ADReal & effective_trial_stress,
                                                     const ADReal & delta_ep)
{
  const ADReal stress_delta =
      effective_trial_stress -
      _elasticity_model->computeMandelStress(delta_ep * _Np[_qp], /*plasticity_update = */ true)
          .doubleContraction(_Np[_qp]);
  const ADReal dstress_delta_ddelta_ep =
      -_elasticity_model->computeMandelStress(_Np[_qp], /*plasticity_update = */ true)
           .doubleContraction(_Np[_qp]);
  const ADReal yield_stress = plasticEnergy(_ep_old[_qp] + delta_ep, 1);
  const ADReal dyield_stress_ddelta_ep = plasticEnergy(_ep_old[_qp] + delta_ep, 2);
  const ADReal dcreep_rate =
      _coefficient * _exponent * std::pow(stress_delta / yield_stress, _exponent - 1);
  return dcreep_rate *
             (dstress_delta_ddelta_ep * yield_stress - stress_delta * dyield_stress_ddelta_ep) /
             yield_stress / yield_stress * _dt -
         1;
}

ADReal
ComputeLargeDeformationXMATStress::plasticEnergy(const ADReal & ep, const unsigned int derivative)
{
  if (derivative == 0)
  {
    _psip_active[_qp] = _sigma_0[_qp] / _arrhenius_coef[_qp] * ep + 0.5 * _eps * ep * ep;
    _psip[_qp] = _gp[_qp] * _psip_active[_qp];
    _dpsip_dd[_qp] = _dgp_dd[_qp] * _psip_active[_qp];
    return _psip[_qp];
  }

  if (derivative == 1)
    return _gp[_qp] * (_sigma_0[_qp] / _arrhenius_coef[_qp] + _eps * ep);

  if (derivative == 2)
    return _gp[_qp] * _eps;

  mooseError(name(), "internal error: unsupported derivative order.");
  return 0;
}

// David's Stuff-----------------------
void
ComputeLargeDeformationXMATStress::computeDelValues(const ADReal & delta_ep)
{
  _del_delt[_qp] = _gp[_qp] * (_sigma_0[_qp] / _arrhenius_coef[_qp]) * delta_ep;
  // return _del_delt[_qp];
}
void
ComputeLargeDeformationXMATStress::initialSetup()
{
  _elasticity_model =
      dynamic_cast<LargeDeformationElasticityModel *>(&getMaterial("elasticity_model"));
  if (!_elasticity_model)
    paramError("elasticity_model",
               "Elasticity model " + _elasticity_model->name() +
                   " is not compatible with ComputeLargeDeformationStress");
}
//---------------------------------
