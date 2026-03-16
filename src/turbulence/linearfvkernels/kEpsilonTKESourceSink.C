//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "kEpsilonTKESourceSink.h"

#include "Assembly.h"
#include "SubProblem.h"
#include "NavierStokesMethods.h"

registerMooseObject("OpenPronghornApp", kEpsilonTKESourceSink);

InputParameters
kEpsilonTKESourceSink::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();

  params.addClassDescription("Elemental kernel to compute the production and destruction "
                             "terms of turbulent kinetic energy (TKE).");

  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");

  params.addRequiredParam<MooseFunctorName>(NS::TKED,
                                            "Coupled turbulent kinetic energy dissipation rate.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Fluid density");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>(NS::mu_t, "Turbulent viscosity.");

  params.addParam<std::vector<BoundaryName>>(
      "walls", {}, "Boundaries that correspond to solid walls.");

  MooseEnum wall_treatment("eq_newton eq_incremental eq_linearized neq", "neq");
  params.addParam<MooseEnum>("wall_treatment", wall_treatment, "Method used for wall functions.");

  params.addParam<Real>("C_mu", 0.09, "Turbulent kinetic energy closure constant C_mu.");
  params.addParam<Real>("C_pl", 10.0, "Production limiter multiplier C_pl.");

  // ---- Extensions for k-epsilon variants and corrections (bulk production) ----

  MooseEnum k_eps_variant("Standard StandardLowRe StandardTwoLayer Realizable RealizableTwoLayer",
                          "Standard");
  params.addParam<MooseEnum>("k_epsilon_variant",
                             k_eps_variant,
                             "k-epsilon model variant used for the bulk production term.");

  params.addParam<bool>(
      "use_buoyancy", false, "Include buoyancy production Gb in the bulk k-equation.");
  params.addParam<bool>("use_compressibility",
                        false,
                        "Include compressibility correction gamma_M in the bulk k-equation.");
  params.addParam<bool>("use_curvature_correction",
                        false,
                        "Apply a curvature correction factor f_c to the shear production.");

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
      "Minimum k value used in the implicit destruction coefficient rho*eps/k to "
      "prevent division by zero when k is near zero during early iterations.");

  params.addParam<Real>(
      "eps_min",
      1e-10,
      "Minimum epsilon used inside the C_pl production limiter "
      "min(Gk, C_pl * rho * max(eps, eps_min)).  Prevents the limiter from "
      "being trivially zero when eps is near zero at initialisation, which would "
      "allow unbounded net TKE growth (C_pl-1 terms accumulate each iteration).");

  params.addParam<Real>(
      "mu_t_prod_max",
      5000.0,
      "Maximum mu_t / mu ratio applied inside the production and wall-shear "
      "computations.  Should match (or be <= to) the mu_t_ratio_max set in "
      "kEpsilonViscosity.  Caps stale over-large mu_t values before k and epsilon "
      "have converged to a consistent state.");

  params.addParam<Real>(
      "eps_functor_max",
      1e6,
      "Maximum epsilon value read from the TKED functor when computing the TKE bulk "
      "destruction coefficient rho*eps/k.  Near-wall cells in the TKED two-layer "
      "region carry large algebraic eps values (physically correct near the wall, "
      "but potentially problematic for immediately adjacent bulk k-cells).  Setting "
      "this to a reasonable physical upper bound (e.g. 1e4–1e6) prevents runaway "
      "destruction in those bulk cells.  Should match or exceed tked_max_phys "
      "used in kEpsilonTKEDSourceSink.");

  params.addParam<Real>(
      "tke_min_phys",
      1e-10,
      "Physical lower bound on TKE (k). When the current cell k falls below this "
      "value a strong penalty source is injected to drive it back to tke_min_phys. "
      "This prevents negative k values from the linear solver causing NaN/Inf in "
      "subsequent iterations.");

  params.addParam<Real>(
      "tke_max_phys",
      1e4,
      "Physical upper bound on TKE (k). When the current cell k exceeds this value "
      "a strong penalty source drives it back. Default 1e4 m²/s² covers most "
      "industrial flows; raise for supersonic/very-high-speed cases.");

  params.addParam<Real>(
      "C_pk",
      0.0,
      "Durbin realizability coefficient for the k-based production limiter: "
      "G_k <= C_pk * rho * k * |S|.  A k-based bound that is effective even "
      "when eps ~ 0 at start-up.  0 = disabled (default).  "
      "Recommended value when enabled: 0.667 (= 2/3).");

  params.addParam<bool>(
      "use_kato_launder",
      false,
      "Replace the standard shear production G_k = mu_t S^2 with the "
      "Kato–Launder (1993) form G_k = mu_t |S| |Omega|. This eliminates the "
      "stagnation-point anomaly and significantly improves stability on coarse "
      "or poorly-conditioned meshes.");

  MooseEnum grad_method_enum("moose_functor local_least_squares", "moose_functor");
  params.addParam<MooseEnum>(
      "gradient_method",
      grad_method_enum,
      "Velocity gradient method used for turbulence production: "
      "'moose_functor' (default, Green-Gauss or whatever the variable uses) or "
      "'local_least_squares' (bypass MOOSE, reconstruct via face-neighbour "
      "inverse-distance-weighted LS — recommended for poor-quality meshes).");

  return params;
}

kEpsilonTKESourceSink::kEpsilonTKESourceSink(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(getFunctor<Real>("u")),
    _v_var(params.isParamValid("v") ? &(getFunctor<Real>("v")) : nullptr),
    _w_var(params.isParamValid("w") ? &(getFunctor<Real>("w")) : nullptr),
    _epsilon(getFunctor<Real>(NS::TKED)),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _mu_t(getFunctor<Real>(NS::mu_t)),
    _wall_boundary_names(getParam<std::vector<BoundaryName>>("walls")),
    _wall_treatment(getParam<MooseEnum>("wall_treatment").getEnum<NS::WallTreatmentEnum>()),
    _C_mu(getParam<Real>("C_mu")),
    _C_pl(getParam<Real>("C_pl")),
    _variant(getParam<MooseEnum>("k_epsilon_variant").getEnum<NS::KEpsilonVariant>()),
    _switches({/*use_buoyancy*/ getParam<bool>("use_buoyancy"),
               /*use_compressibility*/ getParam<bool>("use_compressibility"),
               /*use_yap*/ false,
               /*use_low_re_Gprime*/ false,
               /*use_nonlinear*/ false,
               /*use_curvature_correction*/ getParam<bool>("use_curvature_correction")}),
    _Pr_t(getParam<Real>("Pr_t")),
    _C_M(getParam<Real>("C_M")),
    _g(getParam<RealVectorValue>("gravity")),
    _T_functor(params.isParamValid("temperature") ? &(getFunctor<Real>("temperature")) : nullptr),
    _beta_functor(params.isParamValid("beta") ? &(getFunctor<Real>("beta")) : nullptr),
    _c_functor(params.isParamValid("speed_of_sound") ? &(getFunctor<Real>("speed_of_sound"))
                                                     : nullptr),
    _has_T(params.isParamValid("temperature")),
    _has_beta(params.isParamValid("beta")),
    _has_c(params.isParamValid("speed_of_sound")),
    _nonlinear_model(
        getParam<MooseEnum>("nonlinear_model").getEnum<NS::NonlinearConstitutiveRelation>()),
    _curvature_model(getParam<MooseEnum>("curvature_model").getEnum<NS::CurvatureCorrectionModel>()),
    _k_min(getParam<Real>("k_min")),
    _eps_functor_max(getParam<Real>("eps_functor_max")),
    _eps_min(getParam<Real>("eps_min")),
    _mu_t_prod_max(getParam<Real>("mu_t_prod_max")),
    _C_pk(getParam<Real>("C_pk")),
    _tke_min_phys(getParam<Real>("tke_min_phys")),
    _tke_max_phys(getParam<Real>("tke_max_phys")),
    _use_kato_launder(getParam<bool>("use_kato_launder")),
    _grad_method(getParam<MooseEnum>("gradient_method") == "local_least_squares"
                     ? NS::TurbVelocityGradientMethod::LocalLeastSquares
                     : NS::TurbVelocityGradientMethod::MooseFunctor)
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");

  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied.");

  // Additional parameter consistency / applicability checks

  const bool is_two_layer_variant = (_variant == NS::KEpsilonVariant::StandardTwoLayer ||
                                     _variant == NS::KEpsilonVariant::RealizableTwoLayer);

  auto user_set = [&](const std::string & pname) { return params.isParamSetByUser(pname); };

  // Pull the enums locally (avoid name collisions + ensures they're in-scope)
  const MooseEnum nl_local = getParam<MooseEnum>("nonlinear_model");
  const MooseEnum cm_local = getParam<MooseEnum>("curvature_model");

  // Basic coefficient sanity
  if (_C_mu <= 0.0)
    paramError("C_mu", "C_mu must be > 0.");
  if (_C_pl <= 0.0)
    paramError("C_pl", "C_pl must be > 0.");

  // Two-layer variants require wall information (blending / wall-distance-based behavior)
  if (is_two_layer_variant && _wall_boundary_names.empty())
    paramError("walls",
               "Two-layer k-epsilon variants require a non-empty 'walls' list so wall-bounded "
               "elements and wall distance can be computed.");

  // Curvature correction: only applicable for realizable variants in this kernel
  const bool use_curv_flag = getParam<bool>("use_curvature_correction");
  const bool curv_model_enabled = (cm_local != "none");

  // If user explicitly enables curvature correction but leaves model as none -> error
  if (user_set("use_curvature_correction") && use_curv_flag && !curv_model_enabled)
    paramError("curvature_model",
               "You set use_curvature_correction = true, but curvature_model is 'none'. "
               "Select curvature_model != 'none' (or disable use_curvature_correction).");

  // If user explicitly sets both knobs, require them to agree
  if (user_set("use_curvature_correction") && user_set("curvature_model"))
  {
    if (use_curv_flag != curv_model_enabled)
      paramError("use_curvature_correction",
                 "Inconsistent curvature correction settings: use_curvature_correction must match "
                 "whether curvature_model is 'none' (disabled) or not 'none' (enabled).");
  }

  // Buoyancy production requirements + guard against no-op knobs
  if (_switches.use_buoyancy)
  {
    if (!_has_T)
      paramError("temperature",
                 "use_buoyancy = true requires providing the 'temperature' functor.");
    if (!_has_beta)
      paramError("beta",
                 "use_buoyancy = true requires providing the 'beta' (thermal expansion) functor.");
    if (_Pr_t <= 0.0)
      paramError("Pr_t", "Pr_t must be > 0 when use_buoyancy = true.");
  }
  else
  {
    if (user_set("temperature"))
      paramError("temperature",
                 "Parameter 'temperature' is only applicable when use_buoyancy = true.");
    if (user_set("beta"))
      paramError("beta", "Parameter 'beta' is only applicable when use_buoyancy = true.");
    if (user_set("Pr_t"))
      paramError("Pr_t", "Parameter 'Pr_t' is only applicable when use_buoyancy = true.");
    if (user_set("gravity"))
      paramError("gravity", "Parameter 'gravity' is only applicable when use_buoyancy = true.");
  }

  // Compressibility correction requirements + guard against no-op knobs
  if (_switches.use_compressibility)
  {
    if (!_has_c)
      paramError("speed_of_sound",
                 "use_compressibility = true requires providing the 'speed_of_sound' functor.");
    if (_C_M <= 0.0)
      paramError("C_M", "C_M must be > 0 when use_compressibility = true.");
  }
  else
  {
    if (user_set("speed_of_sound"))
      paramError("speed_of_sound",
                 "Parameter 'speed_of_sound' is only applicable when use_compressibility = true.");
    if (user_set("C_M"))
      paramError("C_M", "Parameter 'C_M' is only applicable when use_compressibility = true.");
  }

  // Error if using cylindrical coordinates - gradients won't be correct
  if (getBlockCoordSystem() == Moose::COORD_RZ)
    mooseError(name(), " is not valid on blocks that use an RZ coordinate system.");

  // Shear-strain-based production requires velocity gradients. Request them
  // when the underlying functor is a cell-based FV variable.
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
kEpsilonTKESourceSink::initialSetup()
{
  LinearFVElementalKernel::initialSetup();

  // Identify wall-bounded elements and populate wall-distance and face data
  NS::getWallBoundedElements(
      _wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _wall_bounded);
  NS::getWallDistance(_wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _dist);
  NS::getElementFaceArgs(_wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _face_infos);
}

Real
kEpsilonTKESourceSink::computeMatrixContribution()
{
  // Matrix part carries the destruction term in the bulk and the fully implicit
  // near-wall production/destruction.

  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());

  // Physical bounds enforcement: if k is out of [tke_min_phys, tke_max_phys] inject
  // a large penalty diagonal so the next iterate is driven back inside the bounds.
  // The matching RHS penalty is in computeRightHandSideContribution.
  const Real k_old = _var.getElemValue(*_current_elem_info, state);
  if (k_old < _tke_min_phys || k_old > _tke_max_phys)
    return _bounds_penalty * _current_elem_volume;

  const Real rho = _rho(elem_arg, state);

  // Near-wall formulation: destruction goes into the matrix (positive diagonal),
  // production goes onto the explicit RHS via computeRightHandSideContribution.
  // Previously both were combined as (destruction - production) in the matrix,
  // which caused a NEGATIVE diagonal whenever wall production > destruction —
  // the primary driver of TKE → ∞ on coarse / poorly-conditioned meshes.
  if (_wall_bounded.find(_current_elem_info->elem()) != _wall_bounded.end())
  {
    const auto [production, destruction] = computeWallTerms(elem_arg, state);
    (void)production; // used in computeRightHandSideContribution
    return destruction * _current_elem_volume;
  }

  // Bulk region: implicit destruction ρ epsilon / k.
  // Cap eps read from the functor by _eps_functor_max so that large algebraic
  // epsilon values in adjacent TKED two-layer cells cannot destroy all TKE in
  // neighboring bulk cells within a single iteration.
  const Real epsilon = std::min(_epsilon(elem_arg, state), _eps_functor_max);
  const Real TKE = _var.getElemValue(*_current_elem_info, state);

  // Guard k > k_min to prevent rho*eps/k from blowing up when k → 0
  const Real destruction_bulk = rho * epsilon / std::max(TKE, _k_min);
  return destruction_bulk * _current_elem_volume;
}

Real
kEpsilonTKESourceSink::computeRightHandSideContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());

  // Physical bounds enforcement (must mirror the check in computeMatrixContribution):
  // Return penalty * target so that the combined system drives k → target.
  {
    const Real k_old = _var.getElemValue(*_current_elem_info, state);
    if (k_old < _tke_min_phys)
      return _bounds_penalty * _tke_min_phys * _current_elem_volume;
    if (k_old > _tke_max_phys)
      return _bounds_penalty * _tke_max_phys * _current_elem_volume;
  }

  // Near-wall cells: production is explicit (on RHS) to keep the matrix diagonal positive.
  // The coefficient-of-k form means RHS = production_coeff * k_current * V.
  if (_wall_bounded.find(_current_elem_info->elem()) != _wall_bounded.end())
  {
    const auto [production_coeff, destruction_coeff] = computeWallTerms(elem_arg, state);
    (void)destruction_coeff; // used in computeMatrixContribution
    const Real TKE = std::max(_var.getElemValue(*_current_elem_info, state), _k_min);
    return production_coeff * TKE * _current_elem_volume;
  }

  // Bulk production term (with optional corrections)
  Real production = computeBulkProduction(elem_arg, state);

  // Production limiter: C_pl * rho * max(eps, eps_min).
  // Using eps_min prevents the limiter from being trivially zero at start-up
  // (eps ≈ 0 at t=0 would allow unbounded net TKE growth otherwise).
  const Real rho = _rho(elem_arg, state);
  const Real eps = _epsilon(elem_arg, state);
  const Real production_limit = _C_pl * rho * std::max(eps, _eps_min);
  production = std::min(production, production_limit);

  return production * _current_elem_volume;
}

Real
kEpsilonTKESourceSink::computeBulkProduction(const Moose::ElemArg & elem_arg,
                                             const Moose::StateArg & state) const
{
  const Real rho = _rho(elem_arg, state);
  const Real mu  = _mu(elem_arg, state);
  // Cap mu_t to _mu_t_prod_max * mu to guard against stale over-large values
  // when k is large but epsilon has not yet caught up (common in early SIMPLE iterations).
  const Real mu_t = std::min(_mu_t(elem_arg, state), _mu_t_prod_max * mu);
  const Real k = _var.getElemValue(*_current_elem_info, state);
  const Real k_safe = std::max(k, _k_min);
  const Real eps = _epsilon(elem_arg, state);

  // Strain/rotation invariants and divergence — use the selected gradient method
  const Elem * elem = elem_arg.elem;
  auto inv = NS::computeStrainRotationInvariantsEx(
      _u_var, _v_var, _w_var, elem, elem_arg, state, _grad_method);

  // Baseline shear production G_k (standard or Kato–Launder form)
  Real Gk = _use_kato_launder
                ? NS::computeGkKatoLaunder(mu_t, inv.S2, inv.W2, rho, k, inv.div_u,
                                           _switches.use_compressibility)
                : NS::computeGk(mu_t, inv.S2, rho, k, inv.div_u,
                                _switches.use_compressibility);

  // Buoyancy production Gb
  Real Gb = 0.0;
  if (_switches.use_buoyancy && _has_T && _has_beta)
  {
    const Real beta = (*_beta_functor)(elem_arg, state);
    const auto grad_T = _T_functor->gradient(elem_arg, state);
    Gb = NS::computeGb(beta, mu_t, _Pr_t, grad_T, _g);
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

  // Curvature / rotation correction factor f_c (Realizable variants)
  Real fc = 1.0;
  if (_curvature_model != NS::CurvatureCorrectionModel::None)
    fc = NS::computeCurvatureFactor(_curvature_model, inv);

  // Compressibility correction gamma_M (subtracted)
  Real gamma_M = 0.0;
  if (_switches.use_compressibility && _has_c)
  {
    const Real c = (*_c_functor)(elem_arg, state);
    gamma_M = NS::computeGammaM(rho, _C_M, k, eps, c);
  }

  // Variant-dependent form of the shear production (fc only affects realizable variants)
  Real shear_prod = 0.0;
  switch (_variant)
  {
    case NS::KEpsilonVariant::Standard:
    case NS::KEpsilonVariant::StandardLowRe:
    case NS::KEpsilonVariant::StandardTwoLayer:
      shear_prod = Gk;
      break;

    case NS::KEpsilonVariant::Realizable:
    case NS::KEpsilonVariant::RealizableTwoLayer:
      shear_prod = fc * Gk;
      break;

    default:
      mooseError("kEpsilonTKESourceSink: unsupported k-epsilon variant.");
  }

  Real total_prod = shear_prod + Gb + Gnl - gamma_M;

  // Durbin (1996) k-based realizability limiter: G_k ≤ C_pk · ρ · k · |S|
  // Effective even when ε ≈ 0 at start-up (unlike the ε-based C_pl limiter).
  // Disabled by default (C_pk = 0); recommended value: 0.667.
  if (_C_pk > 0.0)
  {
    const Real S_mag = std::sqrt(std::max(inv.S2, 0.0));
    const Real durbin_limit = _C_pk * rho * k_safe * S_mag;
    total_prod = std::min(total_prod, durbin_limit);
  }

  return total_prod;
}

std::pair<Real, Real>
kEpsilonTKESourceSink::computeWallTerms(const Moose::ElemArg & elem_arg,
                                        const Moose::StateArg & state) const
{
  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  // Use k_min guard so that sqrt(TKE) and 1/sqrt(TKE) are always finite
  const Real TKE = std::max(_var.getElemValue(*_current_elem_info, state), _k_min);
  // Maximum allowed mu_t for wall-shear computation (prevents stale over-large μ_t values)
  const Real mu_t_cap = _mu_t_prod_max * mu;

  const auto & face_info_vec = libmesh_map_find(_face_infos, _current_elem_info->elem());
  const auto & distance_vec = libmesh_map_find(_dist, _current_elem_info->elem());
  mooseAssert(distance_vec.size(), "Wall-bounded element without distance data.");
  mooseAssert(distance_vec.size() == face_info_vec.size(),
              "Mismatch between face-info and distance vectors.");

  // Cell-centroid velocity
  RealVectorValue velocity(_u_var(elem_arg, state));
  if (_v_var)
    velocity(1) = (*_v_var)(elem_arg, state);
  if (_w_var)
    velocity(2) = (*_w_var)(elem_arg, state);

  const Real tot_weight = static_cast<Real>(distance_vec.size());
  Real production = 0.0;
  Real destruction = 0.0;

  for (unsigned int i = 0; i < distance_vec.size(); ++i)
  {
    const Real d = distance_vec[i];
    mooseAssert(d > 0.0, "Wall distance must be positive.");

    const auto * fi = face_info_vec[i];
    const RealVectorValue tangential_vel =
        velocity - (velocity * fi->normal()) * fi->normal();
    const Real parallel_speed = NS::computeSpeed<Real>(tangential_vel);

    // y+ using either NEQ (non-iterative) or equilibrium (iterative) method
    Real y_plus = 0.0;
    if (_wall_treatment == NS::WallTreatmentEnum::NEQ)
      y_plus = d * std::sqrt(std::sqrt(_C_mu) * TKE) * rho / mu;
    else
      y_plus = NS::findyPlus<Real>(mu, rho, std::max(parallel_speed, 1e-10), d);

    // Evaluate mu_t at the wall face (capped to avoid runaway)
    const bool on_elem_side = _var.hasFaceSide(*fi, true);
    const Elem * loc_elem = on_elem_side ? &fi->elem() : fi->neighborPtr();
    const Moose::FaceArg facearg = {
        fi, Moose::FV::LimiterType::CentralDifference, false, false, loc_elem, nullptr};
    const Real wall_mu_t = std::min(_mu_t(facearg, state), mu_t_cap);
    const Real wall_mu   = _mu(facearg, state);

    // Wall shear stress τ_w = (μ_t + μ) * |du/dy|
    const Real tau_w = (wall_mu_t + wall_mu) * (parallel_speed / d);

    if (y_plus < 11.25)
    {
      // Viscous sublayer: no log-law production; purely viscous destruction coefficient
      //   destruction_coeff = 2 μ / d²  (linearization of ε_wall = 2ν k / d² then ρε/k)
      destruction += 2.0 * wall_mu / Utility::pow<2>(d) / tot_weight;
    }
    else
    {
      // Log-layer destruction coefficient: C_μ^{3/4} ρ √k / (κ d)  — coefficient of k
      const Real k_half = std::sqrt(TKE);
      destruction += std::pow(_C_mu, 0.75) * rho * k_half /
                     (NS::von_karman_constant * d) / tot_weight;

      // Log-layer production coefficient: τ_w C_μ^{1/4} / (√k κ d) — coefficient of k
      // Bounded so it can never exceed the destruction coefficient (physical realizability)
      const Real prod_coeff = tau_w * std::pow(_C_mu, 0.25) /
                              (k_half * NS::von_karman_constant * d) / tot_weight;
      production += prod_coeff;
    }
  }

  return {production, destruction};
}
