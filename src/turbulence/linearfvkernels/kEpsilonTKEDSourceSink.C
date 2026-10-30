//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "kEpsilonTKEDSourceSink.h"

#include "Assembly.h"
#include "SubProblem.h"
#include "NavierStokesMethods.h"
#include "LinearFVTurbulentBlended.h"

registerMooseObject("OpenPronghornApp", kEpsilonTKEDSourceSink);

InputParameters
kEpsilonTKEDSourceSink::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();

  params.addClassDescription("Elemental kernel to compute the production and destruction "
                             "terms of the turbulent dissipation rate epsilon.");

  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");

  params.addRequiredParam<MooseFunctorName>(NS::TKE, "Coupled turbulent kinetic energy (k).");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Fluid density.");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>(NS::mu_t, "Turbulent viscosity.");

  params.addParam<std::vector<BoundaryName>>(
      "walls", {}, "Boundaries that correspond to solid walls.");

  MooseEnum wall_treatment("eq_newton eq_incremental eq_linearized neq", "neq");
  params.addParam<MooseEnum>("wall_treatment", wall_treatment, "Method used for wall functions.");

  params.addParam<Real>("C_mu", 0.09, "Turbulent kinetic energy closure constant C_mu.");

  // ---- Epsilon model constants ----
  params.addParam<Real>("C1_eps", 1.44, "C1 epsilon production coefficient.");
  params.addParam<Real>("C2_eps", 1.92, "C2 epsilon destruction coefficient.");
  params.addParam<Real>(
      "C3_eps", 1.0, "C3 epsilon buoyancy coefficient (may be modified by Gb sign).");
  params.addParam<Real>("Ct", 6.0, "Ct coefficient used in Yap / ambient epsilon source terms.");

  params.addParam<Real>(
      "C_lowRe", 1.0, "Low-Re f2 constant C used in the standard low-Re k-epsilon model.");
  params.addParam<Real>("D_lowRe", 1.0, "Low-Re extra production coefficient D for G'.");
  params.addParam<Real>("E_lowRe", 0.00375, "Low-Re extra production coefficient E for G'.");
  params.addParam<Real>("Cw", 0.83, "Yap correction constant Cw for epsilon equation.");

  // ---- Production limiter (to match old LinearFVTKEDSourceSink) ----
  params.addParam<Real>("C_pl", 10.0, "Production limiter constant multiplier.");

  // ---- k-epsilon variants and switches ----
  MooseEnum k_eps_variant("Standard StandardLowRe StandardTwoLayer Realizable RealizableTwoLayer",
                          "Standard");
  params.addParam<MooseEnum>(
      "k_epsilon_variant", k_eps_variant, "k-epsilon model variant used for epsilon.");

  params.addParam<bool>(
      "use_buoyancy", false, "Include buoyancy production Gb in the epsilon equation.");
  params.addParam<bool>("use_compressibility",
                        false,
                        "Include compressibility correction gamma_M in the epsilon "
                        "formulation (mainly in the k-equation, kept for parity).");
  params.addParam<bool>("use_curvature_correction",
                        false,
                        "Apply a curvature correction factor f_c to the shear production.");
  params.addParam<bool>("use_yap", false, "Include Yap correction gamma_y in epsilon.");
  params.addParam<bool>(
      "use_low_re_Gprime", false, "Include additional low-Re production G' in epsilon.");

  params.addParam<Real>(
      "Pr_t", 0.9, "Turbulent Prandtl number used in the buoyancy production term.");
  params.addParam<Real>("C_M", 1.0, "Coefficient used in the compressibility correction gamma_M.");

  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0.0, 0.0, 0.0), "Gravity vector used for buoyancy production.");

  params.addParam<MooseFunctorName>("temperature",
                                    "Temperature functor used for buoyancy production.");
  params.addParam<MooseFunctorName>(
      "beta", "Thermal expansion coefficient functor used for buoyancy production (beta(T)).");
  params.addParam<MooseFunctorName>(
      "speed_of_sound", "Speed of sound functor used for compressibility correction gamma_M.");
  params.addParam<MooseFunctorName>(
      "wall_distance", "Distance to the closest wall; used for Yap and low-Re corrections.");

  // Two-layer quantities, provided by a kEpsilonTwoLayerFunctorMaterial
  params.addParam<MooseFunctorName>(
      "two_layer_lambda",
      "Wall-proximity blending factor lambda (0 = viscosity affected, 1 = fully turbulent); "
      "required for two-layer variants.");
  params.addParam<MooseFunctorName>(
      "two_layer_epsilon",
      "Algebraic (damped) near-wall epsilon; required for two-layer variants.");

  MooseEnum nonlinear_model("none quadratic cubic", "none");
  params.addParam<MooseEnum>("nonlinear_model",
                             nonlinear_model,
                             "Non-linear constitutive relation for standard k-epsilon.");

  MooseEnum curvature_model("none standard", "none");
  params.addParam<MooseEnum>("curvature_model",
                             curvature_model,
                             "Curvature / rotation correction model: "
                             "'none' (default) or 'standard' (Spalart–Shur style).");

  return params;
}

kEpsilonTKEDSourceSink::kEpsilonTKEDSourceSink(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(getFunctor<Real>("u")),
    _v_var(params.isParamValid("v") ? &(getFunctor<Real>("v")) : nullptr),
    _w_var(params.isParamValid("w") ? &(getFunctor<Real>("w")) : nullptr),
    _k(getFunctor<Real>(NS::TKE)),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _mu_t(getFunctor<Real>(NS::mu_t)),
    _wall_boundary_names(getParam<std::vector<BoundaryName>>("walls")),
    _wall_treatment(getParam<MooseEnum>("wall_treatment").getEnum<NS::WallTreatmentEnum>()),
    _C_mu(getParam<Real>("C_mu")),
    _C1_eps(getParam<Real>("C1_eps")),
    _C2_eps(getParam<Real>("C2_eps")),
    _C3_eps(getParam<Real>("C3_eps")),
    _Ct(getParam<Real>("Ct")),
    _C_lowRe(getParam<Real>("C_lowRe")),
    _D_lowRe(getParam<Real>("D_lowRe")),
    _E_lowRe(getParam<Real>("E_lowRe")),
    _Cw(getParam<Real>("Cw")),
    _C_pl(getParam<Real>("C_pl")),
    _variant(getParam<MooseEnum>("k_epsilon_variant").getEnum<NS::KEpsilonVariant>()),
    _two_layer(_variant == NS::KEpsilonVariant::StandardTwoLayer ||
               _variant == NS::KEpsilonVariant::RealizableTwoLayer),
    _switches{/*use_buoyancy*/ getParam<bool>("use_buoyancy"),
              /*use_compressibility*/ getParam<bool>("use_compressibility"),
              /*use_yap*/ getParam<bool>("use_yap"),
              /*use_low_re_Gprime*/ getParam<bool>("use_low_re_Gprime"),
              /*use_nonlinear*/ false,
              /*use_curvature_correction*/ getParam<bool>("use_curvature_correction")},
    _Pr_t(getParam<Real>("Pr_t")),
    _C_M(getParam<Real>("C_M")),
    _g(getParam<RealVectorValue>("gravity")),
    _T_functor(params.isParamValid("temperature") ? &(getFunctor<Real>("temperature")) : nullptr),
    _beta_functor(params.isParamValid("beta") ? &(getFunctor<Real>("beta")) : nullptr),
    _c_functor(params.isParamValid("speed_of_sound") ? &(getFunctor<Real>("speed_of_sound"))
                                                     : nullptr),
    _wall_distance_functor(
        params.isParamValid("wall_distance") ? &(getFunctor<Real>("wall_distance")) : nullptr),
    _two_layer_lambda(params.isParamValid("two_layer_lambda")
                          ? &(getFunctor<Real>("two_layer_lambda"))
                          : nullptr),
    _two_layer_epsilon(params.isParamValid("two_layer_epsilon")
                           ? &(getFunctor<Real>("two_layer_epsilon"))
                           : nullptr),
    _has_T(params.isParamValid("temperature")),
    _has_beta(params.isParamValid("beta")),
    _has_c(params.isParamValid("speed_of_sound")),
    _has_wall_distance(params.isParamValid("wall_distance")),
    _nonlinear_model(
        getParam<MooseEnum>("nonlinear_model").getEnum<NS::NonlinearConstitutiveRelation>()),
    _curvature_model(getParam<MooseEnum>("curvature_model").getEnum<NS::CurvatureCorrectionModel>())
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");
  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied.");

  // Additional constructor error checks (variant/option consistency)

  const bool is_two_layer = _two_layer;

  const bool is_low_re = (_variant == NS::KEpsilonVariant::StandardLowRe);

  // Helper: "user explicitly set this parameter"
  auto user_set = [&](const std::string & pname) { return params.isParamSetByUser(pname); };

  // Pull enums locally so we don't depend on other names being in-scope
  const MooseEnum nl_local = getParam<MooseEnum>("nonlinear_model");
  const MooseEnum cm_local = getParam<MooseEnum>("curvature_model");

  // ---- Wall-distance requirements (variants + enabled options) ----
  if ((is_two_layer || is_low_re) && !_has_wall_distance)
    paramError("wall_distance",
               "The selected k_epsilon_variant requires 'wall_distance' (two-layer or low-Re).");

  if (is_two_layer)
  {
    if (!_two_layer_lambda)
      paramError("two_layer_lambda",
                 "Two-layer variants require the 'two_layer_lambda' functor "
                 "(see kEpsilonTwoLayerFunctorMaterial).");
    if (!_two_layer_epsilon)
      paramError("two_layer_epsilon",
                 "Two-layer variants require the 'two_layer_epsilon' functor "
                 "(see kEpsilonTwoLayerFunctorMaterial).");
  }
  else
  {
    if (user_set("two_layer_lambda"))
      paramError("two_layer_lambda",
                 "Parameter 'two_layer_lambda' is only applicable to two-layer variants.");
    if (user_set("two_layer_epsilon"))
      paramError("two_layer_epsilon",
                 "Parameter 'two_layer_epsilon' is only applicable to two-layer variants.");
  }

  if (_switches.use_yap && !_has_wall_distance)
    paramError("wall_distance", "'use_yap=true' requires the 'wall_distance' functor.");

  if (_switches.use_low_re_Gprime && !_has_wall_distance)
    paramError("wall_distance", "'use_low_re_Gprime=true' requires the 'wall_distance' functor.");

  // G' only meaningful for: StandardLowRe
  if (_switches.use_low_re_Gprime && !is_low_re)
    paramError("use_low_re_Gprime",
               "'use_low_re_Gprime=true' is only supported for the StandardLowRe variant.");

  // ---- Cross-parameter consistency: use_curvature_correction vs curvature_model ----
  const bool use_curv_flag = getParam<bool>("use_curvature_correction");
  const bool curv_model_enabled = (cm_local != "none");

  if ((user_set("use_curvature_correction") || user_set("curvature_model")) &&
      (use_curv_flag != curv_model_enabled))
  {
    if (use_curv_flag && !curv_model_enabled)
      paramError(
          "curvature_model",
          "use_curvature_correction=true but curvature_model='none'. "
          "Select a curvature_model (e.g. 'standard') or set use_curvature_correction=false.");
    else if (!use_curv_flag && curv_model_enabled)
      paramError("use_curvature_correction",
                 "curvature_model is enabled but use_curvature_correction=false. "
                 "Either set use_curvature_correction=true or set curvature_model='none'.");
  }

  // ---- Dependent-functor checks for enabled physics toggles ----
  if (_switches.use_buoyancy)
  {
    if (!_has_T)
      paramError("temperature", "'use_buoyancy=true' requires the 'temperature' functor.");
    if (!_has_beta)
      paramError("beta", "'use_buoyancy=true' requires the 'beta' functor.");
  }
  else
  {
    if (user_set("temperature"))
      paramError("temperature",
                 "Parameter 'temperature' is only applicable when use_buoyancy = true.");
    if (user_set("beta"))
      paramError("beta", "Parameter 'beta' is only applicable when use_buoyancy = true.");
  }

  if (_switches.use_compressibility)
  {
    if (!_has_c)
      paramError("speed_of_sound",
                 "'use_compressibility=true' requires the 'speed_of_sound' functor.");
  }
  else
  {
    if (user_set("speed_of_sound"))
      paramError("speed_of_sound",
                 "Parameter 'speed_of_sound' is only applicable when use_compressibility = true.");
  }

  // Error if using cylindrical coordinates - gradients won't be correct
  if (getBlockCoordSystem() == Moose::COORD_RZ)
    mooseError(name(), " is not valid on blocks that use an RZ coordinate system.");

  // Strain / rotation invariants require velocity gradients
  if (_u_var.wrapsType<MooseLinearVariableFV<Real>>())
    requestVariableCellGradient(getParam<MooseFunctorName>("u"));
  if (_v_var && _v_var->wrapsType<MooseLinearVariableFV<Real>>())
    requestVariableCellGradient(getParam<MooseFunctorName>("v"));
  if (_w_var && _w_var->wrapsType<MooseLinearVariableFV<Real>>())
    requestVariableCellGradient(getParam<MooseFunctorName>("w"));

  // keep the switches struct up to date
  _switches.use_nonlinear = (_nonlinear_model != NS::NonlinearConstitutiveRelation::None);
  _switches.use_curvature_correction = (_curvature_model != NS::CurvatureCorrectionModel::None);
}

void
kEpsilonTKEDSourceSink::initialSetup()
{
  LinearFVElementalKernel::initialSetup();

  NS::getWallBoundedElements(
      _wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _wall_bounded);
  NS::getWallDistance(_wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _dist);
  NS::getElementFaceArgs(_wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _face_infos);

  if (_two_layer)
  {
    for (const auto & output : {"two_layer_lambda", "two_layer_epsilon"})
      if (!NS::isTwoLayerOutput(_fe_problem, output, getParam<MooseFunctorName>(output), _tid))
        paramError(output, "Must be the '", output, "' of a kEpsilonTwoLayerFunctorMaterial.");
    checkTwoLayerTransportKernels();
  }
  else
    for (const auto kernel : NS::getTransportKernels(_fe_problem, _sys_num, _var_num, _tid))
      if (const auto b = dynamic_cast<const LinearFVTurbulentBlendedInterface *>(kernel);
          b && b->blendingFactorName())
        kernel->paramError(
            "blending_factor", "Only applicable to two-layer variants of '", name(), "'.");
}

void
kEpsilonTKEDSourceSink::checkTwoLayerTransportKernels() const
{
  // In two-layer variants the epsilon row is lambda-blended with the algebraic epsilon, so every
  // transport kernel of epsilon must scale its rows by the same lambda. A plain kernel, or a
  // blended one without a blending factor, keeps the full transport in the viscosity-affected
  // region and silently breaks the blending.
  const auto & lambda_name = getParam<MooseFunctorName>("two_layer_lambda");

  const auto check_kernel = [&](const MooseObject & kernel)
  {
    const auto blended = dynamic_cast<const LinearFVTurbulentBlendedInterface *>(&kernel);
    if (!blended)
      mooseError("Kernel '",
                 kernel.name(),
                 "' (",
                 kernel.type(),
                 ") acts on '",
                 _var.name(),
                 "' but is not a blended turbulence kernel. Two-layer variants require "
                 "'LinearFVTurbulentBlendedAdvection', 'LinearFVTurbulentBlendedDiffusion' and "
                 "'LinearFVTurbulentBlendedTimeDerivative' for epsilon, with "
                 "'blending_factor = ",
                 lambda_name,
                 "'.");
    if (!blended->blendingFactorName())
      kernel.paramError("blending_factor",
                        "Two-layer variants require 'blending_factor' on every epsilon transport "
                        "kernel. Set it to the 'two_layer_lambda' functor of '",
                        name(),
                        "' ('",
                        lambda_name,
                        "').");
    if (*blended->blendingFactorName() != lambda_name)
      kernel.paramError("blending_factor",
                        "The blending factor '",
                        *blended->blendingFactorName(),
                        "' differs from the 'two_layer_lambda' functor '",
                        lambda_name,
                        "' of '",
                        name(),
                        "'. The epsilon transport and source rows will be blended "
                        "inconsistently.");
  };

  for (const auto kernel : NS::getTransportKernels(_fe_problem, _sys_num, _var_num, _tid))
    check_kernel(*kernel);
}

Real
kEpsilonTKEDSourceSink::computeMatrixContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const Elem * elem = _current_elem_info->elem();

  // ---------------------------------------------------------------------------
  // Two-layer: lambda-blend the epsilon equation with the algebraic epsilon in all cells.
  //   row = lambda * (bulk equation) + (1 - lambda) * a_alg * (eps - eps_alg)
  // The transport kernels apply the same lambda to their rows. In wall-bounded cells their
  // weight is 0, so there epsilon is imposed as lambda * eps_wall + (1 - lambda) * eps_alg
  // (unit coefficient, see the classic branch below). lambda -> 1 recovers the wall function
  // of the high-Re variants, lambda -> 0 the algebraic two-layer epsilon.
  // ---------------------------------------------------------------------------
  if (_two_layer)
  {
    if (_wall_bounded.find(elem) != _wall_bounded.end())
      return 1.0;

    const Real lambda = (*_two_layer_lambda)(elem_arg, state);
    return (lambda * computeBulkDestructionCoeff(elem_arg, state) +
            (1.0 - lambda) * computeAlgebraicCoeff(elem_arg, state)) *
           _current_elem_volume;
  }

  // ---------------------------------------------------------------------------
  // Classic near-wall: algebraic epsilon only in wall-bounded cells
  // ---------------------------------------------------------------------------
  // The coefficient is deliberately 1 rather than the cell volume, matching
  // LinearFVTKEDSourceSink. Advection and diffusion skip wall-bounded rows, but the time
  // derivative still adds rho V / dt; a unit coefficient keeps that negligible for small
  // near-wall cells so the wall-function epsilon is effectively imposed. Scaling by V would
  // make both terms O(V) and pin epsilon to its old value instead.
  if (_wall_bounded.find(elem) != _wall_bounded.end())
    return 1.0;

  return computeBulkDestructionCoeff(elem_arg, state) * _current_elem_volume;
}

Real
kEpsilonTKEDSourceSink::computeRightHandSideContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const Elem * elem = _current_elem_info->elem();

  // Two-layer: see computeMatrixContribution()
  if (_two_layer)
  {
    const Real lambda = (*_two_layer_lambda)(elem_arg, state);
    const Real eps_alg = (*_two_layer_epsilon)(elem_arg, state);

    // Wall-bounded: blended wall-function / algebraic epsilon, not volume-scaled
    if (_wall_bounded.find(elem) != _wall_bounded.end())
      return lambda * computeWallEpsilon(elem_arg, state) + (1.0 - lambda) * eps_alg;

    return (lambda * computeBulkProductionRHS(elem_arg, state) +
            (1.0 - lambda) * computeAlgebraicCoeff(elem_arg, state) * eps_alg) *
           _current_elem_volume;
  }

  // Classic near-wall: not volume-scaled, consistent with the unit matrix coefficient
  if (_wall_bounded.find(elem) != _wall_bounded.end())
    return computeWallEpsilon(elem_arg, state);

  return computeBulkProductionRHS(elem_arg, state) * _current_elem_volume;
}

Real
kEpsilonTKEDSourceSink::computeAlgebraicCoeff(const Moose::ElemArg & elem_arg,
                                              const Moose::StateArg & state) const
{
  // a_alg = rho eps_alg / k = rho sqrt(k) / l_eps, i.e. the inverse of the time scale of the
  // algebraic relaxation, which keeps the algebraic row of the same magnitude as the equation
  const Real k = _k(elem_arg, state);
  const Real eps_alg = (*_two_layer_epsilon)(elem_arg, state);
  return _rho(elem_arg, state) * eps_alg / std::max(k, 1e-20);
}

Real
kEpsilonTKEDSourceSink::computeBulkDestructionCoeff(const Moose::ElemArg & elem_arg,
                                                    const Moose::StateArg & state)
{
  // Implicit destruction term C2_eps f2 rho / T_e (per unit volume)
  const Real rho = _rho(elem_arg, state);
  const Real Te = computeTimeScale(elem_arg, state);

  const Real k = _k(elem_arg, state);
  const Real eps = _var.getElemValue(*_current_elem_info, state);
  const Real nu = _mu(elem_arg, state) / rho;

  Real f2 = 1.0;
  if (_variant == NS::KEpsilonVariant::StandardLowRe)
  {
    const Real Re_t = k * k / std::max(nu * eps, 1e-20);
    f2 = NS::f2_SKE_LRe(_C_lowRe, Re_t);
  }
  else if (_variant == NS::KEpsilonVariant::Realizable ||
           _variant == NS::KEpsilonVariant::RealizableTwoLayer)
    f2 = NS::f2_RKE(k, nu, eps);

  return _C2_eps * f2 * rho / std::max(Te, 1e-20);
}

Real
kEpsilonTKEDSourceSink::computeBulkProductionRHS(const Moose::ElemArg & elem_arg,
                                                 const Moose::StateArg & state)
{
  // Production term C1_eps * Pe / T_e (per unit volume)
  const Real Pe = computeBulkPe(elem_arg, state);
  const Real Te = computeTimeScale(elem_arg, state);
  return _C1_eps * Pe / std::max(Te, 1e-20);
}

Real
kEpsilonTKEDSourceSink::computeTimeScale(const Moose::ElemArg & elem_arg,
                                         const Moose::StateArg & state) const
{
  // Classic high-Re time scale T_e = k / eps.
  const Real k = _k(elem_arg, state);
  const Real eps = _var.getElemValue(*_current_elem_info, state);
  return k / std::max(eps, 1e-20);
}

Real
kEpsilonTKEDSourceSink::computeWallEpsilon(const Moose::ElemArg & elem_arg,
                                           const Moose::StateArg & state) const
{
  const auto & elem = *_current_elem_info->elem();
  const auto & distance_vec = libmesh_map_find(_dist, &elem);
  const auto & face_info_vec = libmesh_map_find(_face_infos, &elem);

  mooseAssert(distance_vec.size(), "Wall-bounded element without distance data.");
  mooseAssert(distance_vec.size() == face_info_vec.size(),
              "Mismatch between face-info and distance data size.");

  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  const Real k = _k(elem_arg, state);

  // Velocity at the cell centroid (for equilibrium wall treatments)
  RealVectorValue velocity(_u_var(elem_arg, state), 0.0, 0.0);
  if (_dim >= 2 && _v_var)
    velocity(1) = (*_v_var)(elem_arg, state);
  if (_dim >= 3 && _w_var)
    velocity(2) = (*_w_var)(elem_arg, state);

  Real eps_sum = 0.0;
  Real tot_weight = 0.0;

  for (unsigned int i = 0; i < distance_vec.size(); ++i)
  {
    const Real d = distance_vec[i];
    mooseAssert(d > 0.0, "Wall distance must be positive.");

    Real y_plus = 0.0;

    if (_wall_treatment == NS::WallTreatmentEnum::NEQ)
    {
      // Non-equilibrium / non-iterative: use local k to estimate u_tau
      y_plus = d * std::sqrt(std::sqrt(_C_mu) * k) * rho / mu;
    }
    else
    {
      const auto * fi = face_info_vec[i];
      const RealVectorValue tangential_velocity =
          velocity - (velocity * fi->normal()) * fi->normal();
      const Real parallel_speed = NS::computeSpeed<Real>(tangential_velocity);

      y_plus = NS::findyPlus<Real>(mu, rho, std::max(parallel_speed, 1e-10), d);
    }

    Real eps_i = 0.0;
    if (y_plus < 11.25)
      // viscous sublayer branch from LinearFVTKEDSourceSink
      eps_i = 2.0 * k * mu / rho / (d * d);
    else
      // log-layer epsilon
      eps_i = std::pow(_C_mu, 0.75) * std::pow(k, 1.5) / (NS::von_karman_constant * d);

    eps_sum += eps_i;
    tot_weight += 1.0;
  }

  if (tot_weight <= 0.0)
    return 0.0;

  return eps_sum / tot_weight;
}

Real
kEpsilonTKEDSourceSink::computeBulkPe(const Moose::ElemArg & elem_arg,
                                      const Moose::StateArg & state)
{
  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  const Real mu_t = _mu_t(elem_arg, state);
  const Real k = _k(elem_arg, state);
  const Real eps = _var.getElemValue(*_current_elem_info, state);

  // Invariants of strain, rotation and divergence
  auto inv = NS::computeStrainRotationInvariants(_u_var, _v_var, _w_var, elem_arg, state);

  // Shear production G_k
  Real Gk = NS::computeGk(mu_t,
                          inv.S2,
                          rho,
                          k,
                          inv.div_u,
                          /*include_compressibility_terms*/ _switches.use_compressibility);

  // For Realizable Pe we also need the pure shear production S_k = mu_t S^2
  const Real Sk = mu_t * inv.S2;

  // Buoyancy production G_b
  Real Gb = 0.0;
  Real C3_eps_local = _C3_eps;

  if (_switches.use_buoyancy)
  {
    const Real beta = (*_beta_functor)(elem_arg, state);
    const auto grad_T = _T_functor->gradient(elem_arg, state);

    Gb = NS::computeGb(beta, mu_t, _Pr_t, grad_T, _g);

    // Optional C3_eps adjustment based on buoyancy direction
    if (_g.norm() > 0.0)
    {
      libMesh::VectorValue<Real> velocity(_u_var(elem_arg, state), 0.0, 0.0);
      if (_dim >= 2)
        velocity(1) = (*_v_var)(elem_arg, state);
      if (_dim >= 3)
        velocity(2) = (*_w_var)(elem_arg, state);

      libMesh::VectorValue<Real> g_versor = _g / _g.norm();
      Real v_parallel = velocity * g_versor;
      Real v_perpendicular = (velocity - (v_parallel * g_versor)).norm();
      v_parallel = std::abs(v_parallel);

      if (v_perpendicular > 0.0)
        C3_eps_local = std::tanh(v_parallel / v_perpendicular);
    }
  }

  // Optional extra non-linear production G_nl from pre-coded turbulence methods
  Real Gnl = 0.0;
  if (_nonlinear_model != NS::NonlinearConstitutiveRelation::None)
  {
    auto grad_u = NS::computeVelocityGradient(_u_var, _v_var, _w_var, elem_arg, state);
    Gnl = NS::computeGnl(_nonlinear_model, grad_u, inv, mu_t, k, eps);
  }

  // Curvature / rotation correction factor f_c
  Real fc = 1.0;
  if (_curvature_model != NS::CurvatureCorrectionModel::None)
    fc = NS::computeCurvatureFactor(_curvature_model, inv);

  // Optional low-Re extra production G'
  Real Gprime = 0.0;
  if (_switches.use_low_re_Gprime && _has_wall_distance)
  {
    const Real d = (*_wall_distance_functor)(elem_arg, state);
    const Real nu = mu / rho;
    const Real Re_d = std::sqrt(k) * d / std::max(nu, 1e-12);
    const Real Re_t = k * k / std::max(nu * eps, 1e-20);
    const Real f2_LRe = NS::f2_SKE_LRe(_C_lowRe, Re_t);
    Gprime = NS::computeGprime(_D_lowRe, _E_lowRe, f2_LRe, Gk, mu_t, k, d, Re_d);
  }

  // Yap correction gamma_y => translated into a source term (rho / Ct) * gamma_y
  Real yap_source = 0.0;
  if (_switches.use_yap)
  {
    const Real d = (*_wall_distance_functor)(elem_arg, state);
    const Real nu = mu / rho;

    // Limited turbulence time scale for Yap
    const Real eps_safe = std::max(eps, 1e-20);
    const Real T1 = k / eps_safe;
    const Real T2 = 6.0 * nu / eps_safe;
    const Real T3 = std::pow(T1, 1.625) * T2;
    const Real Teff = std::pow(T3, 1.0 / (1.625 + 1.0));

    // Yap length scale from Teff
    const Real l = std::sqrt(std::max(Teff * eps_safe, 0.0));

    const Real cl = NS::cl_from_Cmu(_C_mu);
    const Real l_eps = cl * d;

    const Real gamma_y = NS::computeGammaY(_Cw, eps, k, l, l_eps);
    yap_source = rho / _Ct * gamma_y;
  }

  // Variant-dependent Pe definition
  Real Pe = 0.0;

  switch (_variant)
  {
    case NS::KEpsilonVariant::Standard:
    {
      // Match old LinearFVTKEDSourceSink for Standard:
      Pe = Gk + Gnl + C3_eps_local * Gb;
      break;
    }

    case NS::KEpsilonVariant::StandardTwoLayer:
      // Pepsilon = Gk + Gnl + C3_eps Gb + rho/Ct gamma_y
      Pe = Gk + Gnl + C3_eps_local * Gb + yap_source;
      break;

    case NS::KEpsilonVariant::StandardLowRe:
      // Pepsilon = Gk + Gnl + G' + C3_eps Gb + rho/Ct gamma_y
      Pe = Gk + Gnl + Gprime + C3_eps_local * Gb + yap_source;
      break;

    case NS::KEpsilonVariant::Realizable:
      // Realizable: Pe = f_c S_k + C3_eps Gb
      Pe = fc * Sk + C3_eps_local * Gb;
      break;

    case NS::KEpsilonVariant::RealizableTwoLayer:
      // Realizable two-layer: same, with optional Yap term
      Pe = fc * Sk + C3_eps_local * Gb + yap_source;
      break;

    default:
      mooseError("kEpsilonTKEDSourceSink: unsupported k-epsilon variant.");
  }

  return Pe;
}
