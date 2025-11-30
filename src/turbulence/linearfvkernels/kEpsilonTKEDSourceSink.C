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

  params.addParam<bool>("linearized_model",
                        true,
                        "If true, this kernel is intended to be used in a linear "
                        "solve (kept for API parity).");

  MooseEnum wall_treatment("eq_newton eq_incremental eq_linearized neq", "neq");
  params.addParam<MooseEnum>("wall_treatment",
                             wall_treatment,
                             "Method used for wall functions: "
                             "'eq_newton', 'eq_incremental', 'eq_linearized', or 'neq'.");

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
  params.addParam<Real>("E_lowRe", 1.0, "Low-Re extra production coefficient E for G'.");
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
  params.addParam<bool>(
      "use_nonlinear", false, "Add an extra non-linear production contribution G_nl.");
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
  params.addParam<MooseFunctorName>("nonlinear_production",
                                    "Optional extra non-linear production term G_nl.");
  params.addParam<MooseFunctorName>("curvature_factor",
                                    "Optional curvature correction factor f_c.");
  params.addParam<MooseFunctorName>(
      "wall_distance", "Distance to the closest wall; used for Yap and low-Re corrections.");

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
    _linearized_model(getParam<bool>("linearized_model")),
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
              /*use_nonlinear*/ getParam<bool>("use_nonlinear"),
              /*use_curvature_correction*/ getParam<bool>("use_curvature_correction")},
    _Pr_t(getParam<Real>("Pr_t")),
    _C_M(getParam<Real>("C_M")),
    _g(getParam<RealVectorValue>("gravity")),
    _T_functor(params.isParamValid("temperature") ? &(getFunctor<Real>("temperature")) : nullptr),
    _beta_functor(params.isParamValid("beta") ? &(getFunctor<Real>("beta")) : nullptr),
    _c_functor(params.isParamValid("speed_of_sound") ? &(getFunctor<Real>("speed_of_sound"))
                                                     : nullptr),
    _Gnl_functor(params.isParamValid("nonlinear_production")
                     ? &(getFunctor<Real>("nonlinear_production"))
                     : nullptr),
    _fc_functor(params.isParamValid("curvature_factor") ? &(getFunctor<Real>("curvature_factor"))
                                                        : nullptr),
    _wall_distance_functor(
        params.isParamValid("wall_distance") ? &(getFunctor<Real>("wall_distance")) : nullptr),
    _has_T(params.isParamValid("temperature")),
    _has_beta(params.isParamValid("beta")),
    _has_c(params.isParamValid("speed_of_sound")),
    _has_Gnl(params.isParamValid("nonlinear_production")),
    _has_fc(params.isParamValid("curvature_factor")),
    _has_wall_distance(params.isParamValid("wall_distance"))
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");
  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied.");

  // Strain / rotation invariants require velocity gradients
  if (dynamic_cast<const MooseLinearVariableFV<Real> *>(&_u_var))
    requestVariableCellGradient(getParam<MooseFunctorName>("u"));
  if (_v_var && dynamic_cast<const MooseLinearVariableFV<Real> *>(_v_var))
    requestVariableCellGradient(getParam<MooseFunctorName>("v"));
  if (_w_var && dynamic_cast<const MooseLinearVariableFV<Real> *>(_w_var))
    requestVariableCellGradient(getParam<MooseFunctorName>("w"));
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

    const Real Red = k * d / std::max(nu, 1e-12);

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
  {
    f2 = NS::f2_RKE(k, nu, eps);
  }

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

    const Real Red = k * d / std::max(nu, 1e-12);

    const Real Rey_star = 60.0; // must match the value used in computeMatrixContribution

    if (Red <= Rey_star)
    {
      // Algebraic two-layer epsilon:
      //   epsilon_2L ≈ C_mu^(3/4) k^(3/2) / (kappa * d)
      // This follows the structure of Eq. (1073): epsilon ∝ k^(3/2)/l_eps.
      const Real eps_2layer =
          std::pow(_C_mu, 0.75) * std::pow(k, 1.5) / (NS::von_karman_constant * std::max(d, 1e-12));

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

  // Shear production G_k (compressibility terms handled in k-equation via gamma_M)
  Real Gk = NS::computeGk(mu_t,
                          inv.S2,
                          rho,
                          k,
                          inv.div_u,
                          /*include_compressibility_terms*/ false);

  // For Realizable Pe we also need the pure shear production S_k = mu_t S^2
  const Real Sk = mu_t * inv.S2;

  // Buoyancy production G_b
  Real Gb = 0.0;
  Real C3_eps_local = _C3_eps;

  if (_switches.use_buoyancy && _has_T && _has_beta)
  {
    const Real beta = (*_beta_functor)(elem_arg, state);
    const auto grad_T_raw = _T_functor->gradient(elem_arg, state);
    libMesh::VectorValue<Real> grad_T(grad_T_raw(0), 0., 0.);
    if (_dim >= 2)
      grad_T(1) = grad_T_raw(1);
    if (_dim >= 3)
      grad_T(2) = grad_T_raw(2);

    Gb = NS::computeGb(beta, mu_t, _Pr_t, grad_T, _g);

    // Optional C3_eps adjustment based on buoyancy direction
    if (_g.norm() > 0.0)
    {
      libMesh::VectorValue<Real> velocity(_u_var(elem_arg, state), 0.0, 0.0);
      if (_dim >= 2 && _v_var)
        velocity(1) = (*_v_var)(elem_arg, state);
      if (_dim >= 3 && _w_var)
        velocity(2) = (*_w_var)(elem_arg, state);

      libMesh::VectorValue<Real> g_versor = _g / _g.norm();
      Real v_parallel = velocity * g_versor;
      Real v_perpendicular = (velocity - (v_parallel * g_versor)).norm();
      v_parallel = std::abs(v_parallel);

      if (v_perpendicular > 0.0)
        C3_eps_local = std::tanh(v_parallel / v_perpendicular);
    }
  }

  // Optional extra non-linear production G_nl (only meaningful for Standard variants)
  Real Gnl = 0.0;
  if (_switches.use_nonlinear && _has_Gnl)
    Gnl = (*_Gnl_functor)(elem_arg, state);

  // Curvature correction factor f_c
  Real fc = 1.0;
  if (_switches.use_curvature_correction && _has_fc)
    fc = (*_fc_functor)(elem_arg, state);

  // Optional low-Re extra production G'
  Real Gprime = 0.0;
  if (_switches.use_low_re_Gprime && _has_wall_distance)
  {
    const Real d = (*_wall_distance_functor)(elem_arg, state);
    const Real nu = mu / rho;
    const Real Re_d = k * d / std::max(nu, 1e-12);
    const Real Re_t = k * k / std::max(nu * eps, 1e-20);
    const Real f2_LRe = NS::f2_SKE_LRe(_C_lowRe, Re_t);
    Gprime = NS::computeGprime(_D_lowRe, _E_lowRe, f2_LRe, Gk, mu_t, k, d, Re_d);
  }

  // Yap correction gamma_y => translated into a source term (rho / Ct) * gamma_y
  Real yap_source = 0.0;
  if (_switches.use_yap && _has_wall_distance)
  {
    const Real d = (*_wall_distance_functor)(elem_arg, state);
    const Real nu = mu / rho;

    // Simple turbulent length scale l and epsilon length scale l_eps
    const Real l = std::pow(k, 1.5) / std::max(eps, 1e-20);
    const Real Re_d = k * d / std::max(nu, 1e-12);
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
      // - Gk = mu_t S²
      // - apply C_pl limiter on Gk (production)
      const Real production_limit = _C_pl * rho * eps;
      const Real Gk_limited = std::min(Gk, production_limit);

      // Legacy Standard model had no Gb/Gnl/etc., but we allow them if enabled.
      Pe = Gk_limited + Gnl + C3_eps_local * Gb;
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
      // Realizable: Pe = f_c S_k + C3_eps Gb  (Eq. 1058)
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
