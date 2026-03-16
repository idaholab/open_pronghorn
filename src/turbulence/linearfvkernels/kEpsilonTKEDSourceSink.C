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

  MooseEnum nonlinear_model("none quadratic cubic", "none");
  params.addParam<MooseEnum>("nonlinear_model",
                             nonlinear_model,
                             "Non-linear constitutive relation for standard k-epsilon.");

  MooseEnum curvature_model("none standard", "none");
  params.addParam<MooseEnum>("curvature_model",
                             curvature_model,
                             "Curvature / rotation correction model: "
                             "'none' (default) or 'standard' (Spalart–Shur style).");

  params.addParam<Real>(
      "k_min",
      1e-8,
      "Minimum k value used to guard the turbulent time scale k/eps.");

  params.addParam<Real>(
      "eps_min",
      1e-10,
      "Minimum epsilon used inside the C_pl production limiter "
      "min(Pe, C_pl * rho * max(eps, eps_min)).  Prevents the limiter from "
      "being trivially zero when epsilon is near zero at initialisation.");

  params.addParam<Real>(
      "mu_t_prod_max",
      1e4,
      "Maximum mu_t / mu ratio applied inside the bulk production computation. "
      "Caps stale over-large mu_t values before k and epsilon have converged.");

  params.addParam<Real>(
      "C_pk",
      0.0,
      "Durbin realizability coefficient for the k-based epsilon production limiter: "
      "Pe <= C_pk * rho * k * |S|.  A k-based bound effective when eps ~ 0 at start-up. "
      "0 = disabled (default).  Recommended when enabled: 0.667 (= 2/3).");

  params.addParam<bool>(
      "use_time_scale_limiter",
      true,
      "Apply a Kolmogorov time scale lower bound T = max(k/eps, C_t*sqrt(nu/eps)). "
      "Prevents the implicit epsilon destruction from becoming unboundedly large in "
      "near-wall or highly-dissipative cells. Strongly recommended for stability.");

  params.addParam<Real>(
      "C_t_kolmogorov",
      6.0,
      "Coefficient for the Kolmogorov lower bound on the turbulent time scale "
      "(used when use_time_scale_limiter = true). Default 6.0 matches STAR-CCM+.");

  params.addParam<bool>(
      "use_kato_launder",
      false,
      "Replace G_k = mu_t S^2 with Kato–Launder G_k = mu_t |S| |Omega|. "
      "Must match the setting used in the kEpsilonTKESourceSink kernel.");

  MooseEnum grad_method_enum2("moose_functor local_least_squares", "moose_functor");
  params.addParam<MooseEnum>(
      "gradient_method",
      grad_method_enum2,
      "Velocity gradient method for turbulence production terms. "
      "Must match the setting used in the kEpsilonTKESourceSink kernel.");

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
    _has_T(params.isParamValid("temperature")),
    _has_beta(params.isParamValid("beta")),
    _has_c(params.isParamValid("speed_of_sound")),
    _has_wall_distance(params.isParamValid("wall_distance")),
    _nonlinear_model(
        getParam<MooseEnum>("nonlinear_model").getEnum<NS::NonlinearConstitutiveRelation>()),
    _curvature_model(getParam<MooseEnum>("curvature_model").getEnum<NS::CurvatureCorrectionModel>()),
    _k_min(getParam<Real>("k_min")),
    _eps_min(getParam<Real>("eps_min")),
    _mu_t_prod_max(getParam<Real>("mu_t_prod_max")),
    _C_pk(getParam<Real>("C_pk")),
    _use_time_scale_limiter(getParam<bool>("use_time_scale_limiter")),
    _C_t_kolmogorov(getParam<Real>("C_t_kolmogorov")),
    _use_kato_launder(getParam<bool>("use_kato_launder")),
    _grad_method(getParam<MooseEnum>("gradient_method") == "local_least_squares"
                     ? NS::TurbVelocityGradientMethod::LocalLeastSquares
                     : NS::TurbVelocityGradientMethod::MooseFunctor)
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");
  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied.");

  // Additional constructor error checks (variant/option consistency)

  const bool is_two_layer = (_variant == NS::KEpsilonVariant::StandardTwoLayer ||
                             _variant == NS::KEpsilonVariant::RealizableTwoLayer);

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
}

Real
kEpsilonTKEDSourceSink::computeMatrixContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const Elem * elem = _current_elem_info->elem();

  // ---------------------------------------------------------------------------
  // Two-layer epsilon: algebraic epsilon in the internal near-wall region
  // ---------------------------------------------------------------------------
  if ((_variant == NS::KEpsilonVariant::StandardTwoLayer ||
       _variant == NS::KEpsilonVariant::RealizableTwoLayer) &&
      _has_wall_distance)
  {
    const Real rho = _rho(elem_arg, state);
    const Real mu = _mu(elem_arg, state);
    const Real nu = mu / rho;
    const Real k = _k(elem_arg, state);
    const Real d = (*_wall_distance_functor)(elem_arg, state);

    const Real Red = std::sqrt(k) * d / std::max(nu, 1e-12);

    // Same wall-distance Reynolds threshold used in the two-layer μ_t blending
    const Real Rey_star = 60.0; // if you change this in kEpsilonViscosity, change it here too

    // Viscosity-dominated two-layer region: impose epsilon algebraically
    if (Red <= Rey_star)
      // Matrix coefficient = 1, but in FV we multiply by the cell volume
      return _current_elem_volume;
  }

  // ---------------------------------------------------------------------------
  // Classic near-wall: algebraic epsilon only in wall-bounded cells
  // (used for non two-layer variants)
  // ---------------------------------------------------------------------------
  if (_wall_bounded.find(elem) != _wall_bounded.end() &&
      !(_variant == NS::KEpsilonVariant::StandardTwoLayer ||
        _variant == NS::KEpsilonVariant::RealizableTwoLayer))
    return _current_elem_volume;

  // ---------------------------------------------------------------------------
  // Bulk: implicit destruction term C2_eps f2 rho / T_e
  // ---------------------------------------------------------------------------
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

  const Real destruction_coeff = _C2_eps * f2 * rho / std::max(Te, 1e-20);
  return destruction_coeff * _current_elem_volume;
}

Real
kEpsilonTKEDSourceSink::computeRightHandSideContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const Elem * elem = _current_elem_info->elem();

  // ---------------------------------------------------------------------------
  // Two-layer epsilon: algebraic epsilon in the internal near-wall region
  // ---------------------------------------------------------------------------
  if ((_variant == NS::KEpsilonVariant::StandardTwoLayer ||
       _variant == NS::KEpsilonVariant::RealizableTwoLayer) &&
      _has_wall_distance)
  {
    const Real rho = _rho(elem_arg, state);
    const Real mu = _mu(elem_arg, state);
    const Real nu = mu / rho;
    const Real k = _k(elem_arg, state);
    const Real d = (*_wall_distance_functor)(elem_arg, state);

    const Real Red = std::sqrt(k) * d / std::max(nu, 1e-12);

    const Real Rey_star = 60.0; // must match the value used in computeMatrixContribution

    if (Red <= Rey_star)
    {
      // Algebraic two-layer epsilon (STAR-CCM+ Eq. 1073):
      //   eps_2L = C_mu^(3/4) * k^(3/2) / l_eps
      //   l_eps  = c_l * d * (1 - exp(-Re_d / (2 * c_l)))
      //   c_l    = kappa * C_mu^{-3/4}  (= 0.42 * C_mu^{-3/4})
      //
      // Using the full l_eps with damping function correctly captures the inner
      // near-wall epsilon profile, which grows larger as d -> 0 (l_eps -> 0 faster
      // than k^(3/2) -> 0). The simplified formula kappa*d underestimates epsilon
      // near the wall because it omits the C_mu^{3/4}/c_l = C_mu^{3/2}/kappa factor.
      const Real c_l = NS::cl_from_Cmu(_C_mu);
      const Real A_eps = 2.0 * c_l;
      const Real l_eps =
          c_l * d * (1.0 - std::exp(-Red / std::max(A_eps, 1e-12)));
      const Real eps_2layer = std::pow(k, 1.5) / std::max(l_eps, 1e-12);

      return eps_2layer * _current_elem_volume;
    }
  }

  // ---------------------------------------------------------------------------
  // Classic near-wall: algebraic epsilon only in wall-bounded cells
  // (used for non two-layer variants)
  // ---------------------------------------------------------------------------
  if (_wall_bounded.find(elem) != _wall_bounded.end() &&
      !(_variant == NS::KEpsilonVariant::StandardTwoLayer ||
        _variant == NS::KEpsilonVariant::RealizableTwoLayer))
  {
    const Real eps_wall = computeWallEpsilon(elem_arg, state);
    return eps_wall * _current_elem_volume;
  }

  // ---------------------------------------------------------------------------
  // Bulk: production term C1_eps * Pe / T_e
  // ---------------------------------------------------------------------------
  const Real Pe = computeBulkPe(elem_arg, state);
  const Real Te = computeTimeScale(elem_arg, state);

  const Real production_eps = _C1_eps * Pe / std::max(Te, 1e-20);
  return production_eps * _current_elem_volume;
}

Real
kEpsilonTKEDSourceSink::computeTimeScale(const Moose::ElemArg & elem_arg,
                                         const Moose::StateArg & state) const
{
  const Real k   = std::max(_k(elem_arg, state), _k_min);
  const Real eps = _var.getElemValue(*_current_elem_info, state);
  const Real nu  = _mu(elem_arg, state) / std::max(_rho(elem_arg, state), 1e-20);

  if (_use_time_scale_limiter)
    return NS::limitTurbTimeScale(k, eps, nu, _C_t_kolmogorov);

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
  const Real mu  = _mu(elem_arg, state);
  // Cap mu_t to avoid stale over-large values before k/eps have converged
  const Real mu_t = std::min(_mu_t(elem_arg, state), _mu_t_prod_max * mu);
  const Real k = _k(elem_arg, state);
  const Real k_safe = std::max(k, _k_min);
  const Real eps = _var.getElemValue(*_current_elem_info, state);

  // Invariants of strain, rotation and divergence — use the selected gradient method
  const Elem * elem = elem_arg.elem;
  auto inv = NS::computeStrainRotationInvariantsEx(
      _u_var, _v_var, _w_var, elem, elem_arg, state, _grad_method);

  // Shear production G_k (standard or Kato–Launder form)
  Real Gk = _use_kato_launder
                ? NS::computeGkKatoLaunder(mu_t, inv.S2, inv.W2, rho, k, inv.div_u,
                                           _switches.use_compressibility)
                : NS::computeGk(mu_t, inv.S2, rho, k, inv.div_u,
                                _switches.use_compressibility);

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
    auto grad_u = (_grad_method == NS::TurbVelocityGradientMethod::LocalLeastSquares)
                      ? NS::computeVelocityGradientLS(_u_var, _v_var, _w_var, elem, state)
                      : NS::computeVelocityGradient(_u_var, _v_var, _w_var, elem_arg, state);
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

  // C_pl production limiter (applied consistently across all variants).
  // Use eps_min so the limiter is never trivially zero at start-up (eps ≈ 0).
  const Real production_limit = _C_pl * rho * std::max(eps, _eps_min);

  // Variant-dependent Pe definition
  Real Pe = 0.0;

  switch (_variant)
  {
    case NS::KEpsilonVariant::Standard:
    {
      const Real Gk_limited = std::min(Gk, production_limit);
      Pe = Gk_limited + Gnl + C3_eps_local * Gb;
      break;
    }

    case NS::KEpsilonVariant::StandardTwoLayer:
    {
      // Apply C_pl limiter to Gk for robustness; Yap term is already bounded.
      const Real Gk_limited = std::min(Gk, production_limit);
      Pe = Gk_limited + Gnl + C3_eps_local * Gb + yap_source;
      break;
    }

    case NS::KEpsilonVariant::StandardLowRe:
    {
      const Real Gk_limited = std::min(Gk, production_limit);
      Pe = Gk_limited + Gnl + Gprime + C3_eps_local * Gb + yap_source;
      break;
    }

    case NS::KEpsilonVariant::Realizable:
    {
      // Realizable: Pe = f_c S_k + C3_eps Gb  (STAR-CCM+ Eq. 1058)
      const Real Sk_limited = std::min(fc * Sk, production_limit);
      Pe = Sk_limited + C3_eps_local * Gb;
      break;
    }

    case NS::KEpsilonVariant::RealizableTwoLayer:
    {
      const Real Sk_limited = std::min(fc * Sk, production_limit);
      Pe = Sk_limited + C3_eps_local * Gb + yap_source;
      break;
    }

    default:
      mooseError("kEpsilonTKEDSourceSink: unsupported k-epsilon variant.");
  }

  // Durbin (1996) k-based realizability limiter applied to the epsilon production:
  //   Pe ≤ C_pk · ρ · k · |S|
  // Effective even when ε ≈ 0 at start-up. Disabled by default (C_pk = 0).
  if (_C_pk > 0.0)
  {
    const Real S_mag = std::sqrt(std::max(inv.S2, 0.0));
    const Real durbin_limit = _C_pk * rho * k_safe * S_mag;
    Pe = std::min(Pe, durbin_limit);
  }

  return Pe;
}
