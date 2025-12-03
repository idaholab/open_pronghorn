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
                             "Non-linear constitutive relation for standard k-epsilon: "
                             "'none', 'quadratic', or 'cubic'.");

  MooseEnum curvature_model("none standard", "none");
  params.addParam<MooseEnum>("curvature_model",
                             curvature_model,
                             "Curvature / rotation correction model: "
                             "'none' (default) or 'standard' (Spalart–Shur style).");

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
    _linearized_model(getParam<bool>("linearized_model")),
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
    _nonlinear_model(NS::NonlinearConstitutiveRelation::None),
    _curvature_model(NS::CurvatureCorrectionModel::None)
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");

  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied.");

  // Shear-strain-based production requires velocity gradients. Request them
  // when the underlying functor is a cell-based FV variable.
  if (dynamic_cast<const MooseLinearVariableFV<Real> *>(&_u_var))
    requestVariableCellGradient(getParam<MooseFunctorName>("u"));
  if (_v_var && dynamic_cast<const MooseLinearVariableFV<Real> *>(_v_var))
    requestVariableCellGradient(getParam<MooseFunctorName>("v"));
  if (_w_var && dynamic_cast<const MooseLinearVariableFV<Real> *>(_w_var))
    requestVariableCellGradient(getParam<MooseFunctorName>("w"));

  // Nonlinear model selection
  const auto nl = getParam<MooseEnum>("nonlinear_model");
  if (nl == "quadratic")
    _nonlinear_model = NS::NonlinearConstitutiveRelation::Quadratic;
  else if (nl == "cubic")
    _nonlinear_model = NS::NonlinearConstitutiveRelation::Cubic;
  else
    _nonlinear_model = NS::NonlinearConstitutiveRelation::None;

  // keep the switches struct up to date
  _switches.use_nonlinear = (_nonlinear_model != NS::NonlinearConstitutiveRelation::None);

  // Curvature model section
  const auto cm = getParam<MooseEnum>("curvature_model");
  if (cm == "standard")
    _curvature_model = NS::CurvatureCorrectionModel::Standard;
  else
    _curvature_model = NS::CurvatureCorrectionModel::None;

  // Keep switches in sync
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
  const Real rho = _rho(elem_arg, state);

  // Near-wall formulation: treat production and destruction implicitly
  if (_wall_bounded.find(_current_elem_info->elem()) != _wall_bounded.end())
  {
    const Real mu = _mu(elem_arg, state);
    const Real TKE = _var.getElemValue(*_current_elem_info, state);

    Real production = 0.0;
    Real destruction = 0.0;
    Real tot_weight = 0.0;

    // Local velocity vector at the cell centroid
    RealVectorValue velocity(_u_var(elem_arg, state));
    if (_v_var)
      velocity(1) = (*_v_var)(elem_arg, state);
    if (_w_var)
      velocity(2) = (*_w_var)(elem_arg, state);

    const auto & face_info_vec = libmesh_map_find(_face_infos, _current_elem_info->elem());
    const auto & distance_vec = libmesh_map_find(_dist, _current_elem_info->elem());
    mooseAssert(distance_vec.size(), "Wall-bounded element without distance data.");
    mooseAssert(distance_vec.size() == face_info_vec.size(),
                "Mismatch between face-info and distance data size.");

    std::vector<Real> y_plus_vec;
    std::vector<Real> vel_grad_norm_vec;
    y_plus_vec.reserve(distance_vec.size());
    vel_grad_norm_vec.reserve(distance_vec.size());

    // Compute y+ and velocity-gradient norm at each wall-adjacent face
    for (unsigned int i = 0; i < distance_vec.size(); ++i)
    {
      const auto * fi = face_info_vec[i];
      const Real d = distance_vec[i];

      const RealVectorValue tangential_velocity =
          velocity - (velocity * fi->normal()) * fi->normal();
      const Real parallel_speed = NS::computeSpeed<Real>(tangential_velocity);

      Real y_plus = 0.0;
      if (_wall_treatment == NS::WallTreatmentEnum::NEQ)
      {
        // Non-equilibrium: u_tau estimated from local TKE
        y_plus = d * std::sqrt(std::sqrt(_C_mu) * TKE) * rho / mu;
      }
      else
      {
        // Equilibrium variants: solve for y+ using standard wall functions
        const Real u_t = std::max(parallel_speed, 1e-10); // avoid divisions by zero in helper
        y_plus = NS::findyPlus<Real>(mu, rho, u_t, d);
      }

      const Real vel_grad_norm = parallel_speed / d;

      y_plus_vec.push_back(y_plus);
      vel_grad_norm_vec.push_back(vel_grad_norm);
      tot_weight += 1.0;
    }

    // Accumulate near-wall production/destruction across all wall faces
    for (unsigned int i = 0; i < y_plus_vec.size(); ++i)
    {
      const Real y_plus = y_plus_vec[i];
      const auto * fi = face_info_vec[i];

      const bool defined_on_elem_side = _var.hasFaceSide(*fi, true);
      const Elem * const loc_elem = defined_on_elem_side ? &fi->elem() : fi->neighborPtr();

      Moose::FaceArg face_arg{
          fi, Moose::FV::LimiterType::CentralDifference, false, false, loc_elem, nullptr};

      const Real wall_mut = _mu_t(face_arg, state);
      const Real wall_mu = _mu(face_arg, state);

      const Real tau_w = (wall_mut + wall_mu) * vel_grad_norm_vec[i];
      const Real d = distance_vec[i];

      const Real destruction_visc = 2.0 * wall_mu / (d * d) / tot_weight;
      const Real destruction_log =
          std::pow(_C_mu, 0.75) * rho * std::sqrt(TKE) / (NS::von_karman_constant * d) / tot_weight;

      if (y_plus < 11.25)
        destruction += destruction_visc;
      else
      {
        destruction += destruction_log;
        production += tau_w * std::pow(_C_mu, 0.25) / std::sqrt(TKE) /
                      (NS::von_karman_constant * d) / tot_weight;
      }
    }

    // Implicit near-wall production/destruction (multiplied by TKE when applied)
    return (destruction - production) * _current_elem_volume;
  }

  // Bulk region: implicit destruction ρ epsilon / k
  const Real epsilon = _epsilon(elem_arg, state);
  const Real TKE = _var.getElemValue(*_current_elem_info, state);

  const Real destruction_bulk = rho * epsilon / TKE;
  return destruction_bulk * _current_elem_volume;
}

Real
kEpsilonTKESourceSink::computeRightHandSideContribution()
{
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());

  // No explicit RHS term in near-wall cells: everything handled in the matrix
  if (_wall_bounded.find(_current_elem_info->elem()) != _wall_bounded.end())
    return 0.0;

  // Bulk production term (with optional corrections)
  Real production = computeBulkProduction(elem_arg, state);

  // Production limiter to avoid excessive TKE in stagnation zones
  const Real rho = _rho(elem_arg, state);
  const Real eps = _epsilon(elem_arg, state);
  const Real production_limit = _C_pl * rho * eps;
  production = std::min(production, production_limit);

  return production * _current_elem_volume;
}

Real
kEpsilonTKESourceSink::computeBulkProduction(const Moose::ElemArg & elem_arg,
                                             const Moose::StateArg & state) const
{
  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  const Real mu_t = _mu_t(elem_arg, state);
  const Real k = _var.getElemValue(*_current_elem_info, state);
  const Real eps = _epsilon(elem_arg, state);

  // Strain/rotation invariants and divergence
  auto inv = NS::computeStrainRotationInvariants(_u_var, _v_var, _w_var, elem_arg, state);

  // Baseline shear production G_k
  Real Gk = NS::computeGk(mu_t, inv.S2, rho, k, inv.div_u, /*include_compressibility_terms*/ false);

  // Buoyancy production Gb
  Real Gb = 0.0;
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
  }

  // Optional extra non-linear production G_nl from pre-coded turbulence methods
  Real Gnl = 0.0;
  if (_nonlinear_model != NS::NonlinearConstitutiveRelation::None)
  {
    auto grad_u = NS::computeVelocityGradient(_u_var, _v_var, _w_var, elem_arg, state);
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

  // Total bulk production passed to the RHS
  return shear_prod + Gb + Gnl - gamma_M;
}
