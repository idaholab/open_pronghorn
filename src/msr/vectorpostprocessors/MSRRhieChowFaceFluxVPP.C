#include "MSRRhieChowFaceFluxVPP.h"

#include "RhieChowFaceFluxProvider.h"
#include "MooseMesh.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <utility>
#include <vector>

registerMooseObject("OpenPronghornApp", MSRRhieChowFaceFluxVPP);

InputParameters
MSRRhieChowFaceFluxVPP::validParams()
{
  auto params = GeneralVectorPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "Rhie-Chow face-flux provider to sample.");
  params.addParam<bool>(
      "internal_only",
      true,
      "If true, record only internal faces. This is the intended bridge-validation mode.");
  params.addClassDescription(
      "Records Rhie-Chow volumetric face fluxes for direct linear-vs-reconstructed validation.");
  return params;
}

MSRRhieChowFaceFluxVPP::MSRRhieChowFaceFluxVPP(const InputParameters & params)
  : GeneralVectorPostprocessor(params),
    _rc(getUserObject<RhieChowFaceFluxProvider>("rhie_chow_user_object")),
    _internal_only(getParam<bool>("internal_only")),
    _face_id(declareVector("face_id")),
    _flux(declareVector("volumetric_flux"))
{
}

void
MSRRhieChowFaceFluxVPP::initialize()
{
  _face_id.clear();
  _flux.clear();
}

void
MSRRhieChowFaceFluxVPP::execute()
{
  // This diagnostic is deliberately used only in the small serial smoke test.
  // In serial, face IDs are directly comparable between the identical parent/child meshes.
  std::vector<std::pair<dof_id_type, Real>> samples;
  samples.reserve(_fe_problem.mesh().faceInfo().size());

  for (const auto * fi : _fe_problem.mesh().faceInfo())
  {
    if (!fi)
      continue;

    if (_internal_only && !fi->neighborPtr())
      continue;

    const Real q = _rc.getVolumetricFaceFlux(
        Moose::FV::InterpMethod::RhieChow, *fi, Moose::currentState(), 0, false);

    samples.emplace_back(fi->id(), q);
  }

  std::sort(samples.begin(), samples.end(), [](const auto & a, const auto & b) {
    return a.first < b.first;
  });

  Real sum = 0.0;
  Real sum_abs = 0.0;
  Real sum_sq = 0.0;
  Real weighted_id_sum = 0.0;
  Real max_abs = 0.0;

  for (const auto & sample : samples)
  {
    const auto id = sample.first;
    const Real q = sample.second;

    _face_id.push_back(static_cast<Real>(id));
    _flux.push_back(q);

    sum += q;
    sum_abs += std::abs(q);
    sum_sq += q * q;
    weighted_id_sum += static_cast<Real>(id + 1) * q;
    max_abs = std::max(max_abs, std::abs(q));
  }

  _console << std::setprecision(17)
           << "MSR_RC_FACE_FLUX_SUMMARY name=" << name()
           << " time=" << _t
           << " count=" << samples.size()
           << " sum=" << sum
           << " sum_abs=" << sum_abs
           << " sum_sq=" << sum_sq
           << " weighted_id_sum=" << weighted_id_sum
           << " max_abs=" << max_abs
           << std::endl;
}
