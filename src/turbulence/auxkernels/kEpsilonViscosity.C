//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "kEpsilonViscosity.h"
#include "NavierStokesMethods.h"
#include "NonlinearSystemBase.h"
#include "libmesh/nonlinear_solver.h"

registerMooseObject("OpenPronghornApp", kEpsilonViscosity);

InputParameters
kEpsilonViscosity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Calculates the turbulent viscosity according to the k-epsilon family of models.");

  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::TKE, "Coupled turbulent kinetic energy.");
  params.addRequiredParam<MooseFunctorName>(NS::TKED,
                                            "Coupled turbulent kinetic energy dissipation rate.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity.");

  params.addParam<Real>("C_mu", 0.09, "Base turbulent kinetic energy closure constant C_mu.");
  params.addParam<Real>("mu_t_ratio_max", 1e5, "Maximum allowable mu_t_ratio : mu/mu_t.");

  params.addParam<std::vector<BoundaryName>>(
      "walls", {}, "Boundaries that correspond to solid walls.");

  params.addParam<bool>(
      "bulk_wall_treatment", false, "If true, use classical wall functions for near-wall cells.");

  MooseEnum wall_treatment("eq_newton eq_incremental eq_linearized neq", "eq_newton");
  params.addParam<MooseEnum>(
      "wall_treatment", wall_treatment, "The method used for computing the wall functions.");

  MooseEnum scale_limiter("none standard", "standard");
  params.addParam<MooseEnum>("scale_limiter",
                             scale_limiter,
                             "The method used to limit the k-epsilon time scale: "
                             "'none' or 'standard' (max(T_e, sqrt(nu/epsilon))).");

  // ---- Extended k-epsilon model options ----

  MooseEnum k_eps_variant("Standard StandardLowRe StandardTwoLayer Realizable RealizableTwoLayer",
                          "Standard");
  params.addParam<MooseEnum>(
      "k_epsilon_variant", k_eps_variant, "k-epsilon model variant used for viscosity.");

  MooseEnum two_layer_flavor("Wolfstein NorrisReynolds Xu", "Wolfstein");
  params.addParam<MooseEnum>(
      "two_layer_flavor",
      two_layer_flavor,
      "Two-layer formulation to use for 2L variants (Wolfstein, NorrisReynolds, Xu).");

  // Low-Re damping coefficients for f_mu (SKE LRe)
  params.addParam<Real>("Cd0",
                        0.091,
                        "Low-Re coefficient Cd0 used in the f_mu damping function "
                        "for the Standard k-epsilon Low-Re model.");
  params.addParam<Real>("Cd1",
                        0.0042,
                        "Low-Re coefficient Cd1 used in the f_mu damping function "
                        "for the Standard k-epsilon Low-Re model.");
  params.addParam<Real>("Cd2",
                        0.00011,
                        "Low-Re coefficient Cd2 used in the f_mu damping function "
                        "for the Standard k-epsilon Low-Re model.");

  // Realizable C_mu coefficients (Ca0..Ca3)
  params.addParam<Real>(
      "Ca0", 0.667, "Realizable k-epsilon coefficient Ca0 used in the C_mu expression.");
  params.addParam<Real>(
      "Ca1", 1.25, "Realizable k-epsilon coefficient Ca1 used in the C_mu expression.");
  params.addParam<Real>(
      "Ca2", 1.0, "Realizable k-epsilon coefficient Ca2 used in the C_mu expression.");
  params.addParam<Real>(
      "Ca3", 0.9, "Realizable k-epsilon coefficient Ca3 used in the C_mu expression.");

  // Time-scale coefficient C_t for T = max(Te, C_t sqrt(nu/epsilon))
  params.addParam<Real>("Ct",
                        1.0,
                        "Time-scale constant Ct used in the k-epsilon time scale "
                        "T = max(Te, Ct*sqrt(nu/epsilon)) when scale_limiter='standard'.");

  // Wall distance functor for two-layer and low-Re models
  params.addParam<MooseFunctorName>(
      "wall_distance",
      "Distance to the closest wall; required for two-layer and Low-Re k-epsilon variants.");

  return params;
}

kEpsilonViscosity::kEpsilonViscosity(const InputParameters & params)
  : AuxKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(getFunctor<Real>("u")),
    _v_var(params.isParamValid("v") ? &(getFunctor<Real>("v")) : nullptr),
    _w_var(params.isParamValid("w") ? &(getFunctor<Real>("w")) : nullptr),
    _k(getFunctor<Real>(NS::TKE)),
    _epsilon(getFunctor<Real>(NS::TKED)),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _C_mu(getParam<Real>("C_mu")),
    _mu_t_ratio_max(getParam<Real>("mu_t_ratio_max")),
    _wall_boundary_names(getParam<std::vector<BoundaryName>>("walls")),
    _bulk_wall_treatment(getParam<bool>("bulk_wall_treatment")),
    _wall_treatment(getParam<MooseEnum>("wall_treatment").getEnum<NS::WallTreatmentEnum>()),
    _scale_limiter(getParam<MooseEnum>("scale_limiter")),
    _variant(getParam<MooseEnum>("k_epsilon_variant").getEnum<NS::KEpsilonVariant>()),
    _two_layer_flavor(getParam<MooseEnum>("two_layer_flavor").getEnum<NS::TwoLayerFlavor>()),
    _Cd0(getParam<Real>("Cd0")),
    _Cd1(getParam<Real>("Cd1")),
    _Cd2(getParam<Real>("Cd2")),
    _Ca0(getParam<Real>("Ca0")),
    _Ca1(getParam<Real>("Ca1")),
    _Ca2(getParam<Real>("Ca2")),
    _Ca3(getParam<Real>("Ca3")),
    _Ct(getParam<Real>("Ct")),
    _wall_distance_functor(
        params.isParamValid("wall_distance") ? &(getFunctor<Real>("wall_distance")) : nullptr),
    _has_wall_distance(params.isParamValid("wall_distance"))
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");
  if (_dim >= 3 && !_w_var)
    paramError("w", "In three dimensions, the w velocity must be supplied.");

  // Ensure wall distance is available for variants that need it
  if ((_variant == NS::KEpsilonVariant::StandardTwoLayer ||
       _variant == NS::KEpsilonVariant::RealizableTwoLayer ||
       _variant == NS::KEpsilonVariant::StandardLowRe) &&
      !_has_wall_distance)
    paramError("wall_distance", "The selected k-epsilon variant requires a wall_distance functor.");

  // Variant-specific parameter consistency checks

  const bool is_two_layer_variant = (_variant == NS::KEpsilonVariant::StandardTwoLayer ||
                                     _variant == NS::KEpsilonVariant::RealizableTwoLayer);

  const bool is_realizable_variant = (_variant == NS::KEpsilonVariant::Realizable ||
                                      _variant == NS::KEpsilonVariant::RealizableTwoLayer);

  const bool is_lowre_variant = (_variant == NS::KEpsilonVariant::StandardLowRe);

  // Helper: if user explicitly set a parameter that is not applicable, throw
  auto check_param_applicability =
      [&](const std::string & pname, const bool allowed, const std::string & msg)
  {
    if (params.isParamSetByUser(pname) && !allowed)
      paramError(pname, msg);
  };

  // ---- Two-layer flavor only makes sense for two-layer variants ----
  check_param_applicability(
      "two_layer_flavor",
      is_two_layer_variant,
      "Parameter 'two_layer_flavor' is only applicable when 'k_epsilon_variant' is "
      "'StandardTwoLayer' or 'RealizableTwoLayer'.");

  // ---- Low-Re damping coefficients only apply to StandardLowRe ----
  check_param_applicability(
      "Cd0",
      is_lowre_variant,
      "Parameter 'Cd0' is only applicable when 'k_epsilon_variant' is 'StandardLowRe'.");
  check_param_applicability(
      "Cd1",
      is_lowre_variant,
      "Parameter 'Cd1' is only applicable when 'k_epsilon_variant' is 'StandardLowRe'.");
  check_param_applicability(
      "Cd2",
      is_lowre_variant,
      "Parameter 'Cd2' is only applicable when 'k_epsilon_variant' is 'StandardLowRe'.");

  // ---- Realizable coefficients only apply to realizable variants ----
  check_param_applicability("Ca0",
                            is_realizable_variant,
                            "Parameter 'Ca0' is only applicable when 'k_epsilon_variant' is "
                            "'Realizable' or 'RealizableTwoLayer'.");
  check_param_applicability("Ca1",
                            is_realizable_variant,
                            "Parameter 'Ca1' is only applicable when 'k_epsilon_variant' is "
                            "'Realizable' or 'RealizableTwoLayer'.");
  check_param_applicability("Ca2",
                            is_realizable_variant,
                            "Parameter 'Ca2' is only applicable when 'k_epsilon_variant' is "
                            "'Realizable' or 'RealizableTwoLayer'.");
  check_param_applicability("Ca3",
                            is_realizable_variant,
                            "Parameter 'Ca3' is only applicable when 'k_epsilon_variant' is "
                            "'Realizable' or 'RealizableTwoLayer'.");

  // Error if using cylindrical coordinates - gradients won't be correct
  if (getBlockCoordSystem() == Moose::COORD_RZ)
    mooseError(name(), " is not valid on blocks that use an RZ coordinate system.");

  // Realizable variant needs velocity gradients for strain/rotation invariants
  if (_variant == NS::KEpsilonVariant::Realizable ||
      _variant == NS::KEpsilonVariant::RealizableTwoLayer)
  {
    if (_u_var.wrapsType<MooseLinearVariableFV<Real>>())
      dynamic_cast<MooseLinearVariableFV<Real> *>(&_subproblem.getStandardVariable(_tid, "u"))
          ->computeCellGradients();
    if (_v_var && _v_var->wrapsType<MooseLinearVariableFV<Real>>())
      dynamic_cast<MooseLinearVariableFV<Real> *>(&_subproblem.getStandardVariable(_tid, "v"))
          ->computeCellGradients();
    if (_w_var && _w_var->wrapsType<MooseLinearVariableFV<Real>>())
      dynamic_cast<MooseLinearVariableFV<Real> *>(&_subproblem.getStandardVariable(_tid, "w"))
          ->computeCellGradients();
  }
}

void
kEpsilonViscosity::initialSetup()
{
  if (!_wall_boundary_names.empty())
  {
    NS::getWallBoundedElements(
        _wall_boundary_names, _c_fe_problem, _subproblem, blockIDs(), _wall_bounded);
    NS::getWallDistance(_wall_boundary_names, _c_fe_problem, _subproblem, blockIDs(), _dist);
    NS::getElementFaceArgs(
        _wall_boundary_names, _c_fe_problem, _subproblem, blockIDs(), _face_infos);
  }
}

Real
kEpsilonViscosity::computeValue()
{
  const Elem & elem = *_current_elem;
  const auto elem_arg = makeElemArg(_current_elem);
  const Moose::StateArg state = determineState();

  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  const Real nu = mu / rho;

  // Determine if the element is wall bounded (for bulk wall functions)
  const bool wall_bounded = _wall_bounded.find(_current_elem) != _wall_bounded.end();

  Real mu_t = 0.0;

  // ---------------------------------------------------------------------------
  // Near-wall treatment using classic wall functions (bulk_wall_treatment=true)
  // ---------------------------------------------------------------------------
  if (wall_bounded && _bulk_wall_treatment)
  {
    const auto & elem_distances = _dist.at(&elem);
    const auto min_it = std::min_element(elem_distances.begin(), elem_distances.end());
    const Real min_wall_dist = *min_it;
    const size_t min_index = std::distance(elem_distances.begin(), min_it);
    const auto loc_normal = _face_infos.at(&elem)[min_index]->normal();

    // Local velocity vector
    RealVectorValue velocity(_u_var(elem_arg, state));
    if (_v_var)
      velocity(1) = (*_v_var)(elem_arg, state);
    if (_w_var)
      velocity(2) = (*_w_var)(elem_arg, state);

    // Speed parallel to the wall
    const Real parallel_speed =
        NS::computeSpeed<Real>(velocity - velocity * loc_normal * loc_normal);

    Real y_plus = 0.0;
    Real mut_log = 0.0;
    Real mu_wall = mu;

    if (_wall_treatment == NS::WallTreatmentEnum::EQ_NEWTON)
    {
      const Real u_tau = NS::findUStar<Real>(mu, rho, parallel_speed, min_wall_dist);
      y_plus = min_wall_dist * u_tau * rho / mu;
      mu_wall = rho * Utility::pow<2>(u_tau) * min_wall_dist / parallel_speed;
      mut_log = mu_wall - mu;
    }
    else if (_wall_treatment == NS::WallTreatmentEnum::EQ_INCREMENTAL)
    {
      y_plus = NS::findyPlus<Real>(mu, rho, std::max(parallel_speed, 1e-10), min_wall_dist);
      mu_wall = mu * (NS::von_karman_constant * y_plus /
                      std::log(std::max(NS::E_turb_constant * y_plus, 1 + 1e-4)));
      mut_log = mu_wall - mu;
    }
    else if (_wall_treatment == NS::WallTreatmentEnum::EQ_LINEARIZED)
    {
      const Real a_c = 1.0 / NS::von_karman_constant;
      const Real b_c = 1.0 / NS::von_karman_constant *
                       (std::log(NS::E_turb_constant * std::max(min_wall_dist, 1.0) / mu) + 1.0);
      const Real c_c = parallel_speed;

      const Real u_tau = (-b_c + std::sqrt(b_c * b_c + 4.0 * a_c * c_c)) / (2.0 * a_c);
      y_plus = min_wall_dist * u_tau * rho / mu;
      mu_wall = rho * Utility::pow<2>(u_tau) * min_wall_dist / parallel_speed;
      mut_log = mu_wall - mu;
    }
    else if (_wall_treatment == NS::WallTreatmentEnum::NEQ)
    {
      y_plus = min_wall_dist * std::sqrt(std::sqrt(_C_mu) * _k(elem_arg, state)) * rho / mu;
      mu_wall = mu * (NS::von_karman_constant * y_plus /
                      std::log(std::max(NS::E_turb_constant * y_plus, 1 + 1e-4)));
      mut_log = mu_wall - mu;
    }
    else
      mooseError("kEpsilonViscosity: invalid wall_treatment for bulk wall treatment.");

    if (y_plus <= 5.0)
      mu_t = 0.0; // viscous sublayer
    else if (y_plus >= 30.0)
      mu_t = std::max(mut_log, NS::mu_t_low_limit); // log layer
    else
    {
      // buffer layer blending between viscous and log law
      const Real blending_function = (y_plus - 5.0) / 25.0;
      const Real mut_log_30 = mu * (NS::von_karman_constant * 30.0 /
                                        std::log(std::max(NS::E_turb_constant * 30.0, 1 + 1e-4)) -
                                    1.0);
      mu_t = std::max(blending_function * mut_log_30, NS::mu_t_low_limit);
    }
  }
  // ---------------------------------------------------------------------------
  // Bulk turbulence viscosity from k-epsilon variants
  // ---------------------------------------------------------------------------
  else
  {
    const Real k = _k(elem_arg, state);
    const Real eps = _epsilon(elem_arg, state);

    // Large-eddy time scale T_e = k / epsilon
    const Real Te = k / std::max(eps, 1e-20);

    // Time scale with optional limiter:
    //   T = max(Te, Ct * sqrt(nu/eps))  (STAR-CCM+ Eq. 1046)
    Real time_scale = Te;
    if (_scale_limiter == "standard")
    {
      const Real Tr = _Ct * std::sqrt(nu / std::max(eps, 1e-20));
      time_scale = std::max(Te, Tr);
    }

    // Effective C_mu including damping / realizable modifications
    Real Cmu_eff = _C_mu;

    if (_variant == NS::KEpsilonVariant::StandardLowRe)
    {
      const Real d = (*_wall_distance_functor)(elem_arg, state);
      const Real Red = std::sqrt(k) * d / std::max(nu, 1e-12);
      const Real fmu = NS::fmu_SKE_LRe(_Cd0, _Cd1, _Cd2, Red);
      Cmu_eff = _C_mu * fmu;
    }
    else if (_variant == NS::KEpsilonVariant::Realizable ||
             _variant == NS::KEpsilonVariant::RealizableTwoLayer)
    {
      auto inv = NS::computeStrainRotationInvariants(_u_var, _v_var, _w_var, elem_arg, state);
      Cmu_eff = NS::Cmu_realizable(_Ca0, _Ca1, _Ca2, _Ca3, inv.S2, inv.W2, k, std::max(eps, 1e-20));
    }

    // Base k-epsilon turbulent viscosity
    const Real mu_t_ke = rho * Cmu_eff * k * time_scale;

    // Two-layer blending if requested
    if (_variant == NS::KEpsilonVariant::StandardTwoLayer ||
        _variant == NS::KEpsilonVariant::RealizableTwoLayer)
    {
      mooseAssert(_has_wall_distance, "Two-layer variants require wall_distance functor.");
      const Real d = (*_wall_distance_functor)(elem_arg, state);
      const Real Red = std::sqrt(k) * d / std::max(nu, 1e-12);

      // Two-layer mu ratio and epsilon length scale
      NS::TwoLayerLengths tl;
      if (_two_layer_flavor == NS::TwoLayerFlavor::Wolfstein)
        tl = NS::twoLayerWolfstein(_C_mu, d, Red);
      else if (_two_layer_flavor == NS::TwoLayerFlavor::NorrisReynolds)
        tl = NS::twoLayerNorrisReynolds(_C_mu, d, Red);
      else // Xu (natural convection)
      {
        const Real yv_star = Red * nu / std::max(k, 1e-20);
        tl = NS::twoLayerXu(_C_mu, d, Red, yv_star);
      }

      const Real mu_2layer = tl.mu_ratio * mu;

      // Wall-proximity indicator lambda (0..1)
      const Real Rey_star = 60.0;
      const Real DeltaRey = 10.0;
      const Real A = std::abs(DeltaRey / std::atanh(0.98));
      const Real lambda = 0.5 * (1.0 + std::tanh((Red - Rey_star) / A));

      mu_t = lambda * mu_t_ke + (1.0 - lambda) * mu_2layer;
    }
    else
      mu_t = mu_t_ke;

    mu_t = std::max(mu_t, NS::mu_t_low_limit);
  }

  // Turbulent viscosity limiter
  return std::min(mu_t, _mu_t_ratio_max * mu);
}
