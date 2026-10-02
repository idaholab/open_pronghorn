#pragma once

#include "GeneralVectorPostprocessor.h"

class RhieChowFaceFluxProvider;

/**
 * Diagnostic VectorPostprocessor for comparing Rhie-Chow volumetric face fluxes.
 *
 * Intended for the small serial bridge-validation cases. It records all internal
 * face IDs and volumetric face fluxes, and prints a compact deterministic summary.
 *
 * The input UserObject may be either the linear SIMPLE RhieChowMassFlux family or
 * the nonlinear INSFVRhieChowInterpolator family because both implement
 * RhieChowFaceFluxProvider.
 */
class MSRRhieChowFaceFluxVPP : public GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();
  MSRRhieChowFaceFluxVPP(const InputParameters & params);

  void initialize() override;
  void execute() override;

protected:
  const RhieChowFaceFluxProvider & _rc;
  const bool _internal_only;

  VectorPostprocessorValue & _face_id;
  VectorPostprocessorValue & _flux;
};
