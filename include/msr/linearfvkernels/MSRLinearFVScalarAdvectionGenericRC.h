#pragma once

#include "LinearFVFluxKernel.h"
#include "LinearFVAdvectionDiffusionBC.h"
#include "RhieChowFaceFluxProvider.h"

/**
 * Linear FV passive-scalar advection using the generic RhieChowFaceFluxProvider API.
 *
 * This mirrors MOOSE LinearFVScalarAdvection, but accepts either:
 *   - linear SIMPLE RhieChowMassFlux-derived providers, or
 *   - nonlinear INSFVRhieChowInterpolator-derived providers.
 *
 * That lets a linear BDF2 scalar system use the same reconstructed Rhie-Chow
 * face flux validated for the EPI3V transport child.
 */
class MSRLinearFVScalarAdvectionGenericRC : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();
  MSRLinearFVScalarAdvectionGenericRC(const InputParameters & params);

  Real computeElemMatrixContribution() override;
  Real computeNeighborMatrixContribution() override;
  Real computeElemRightHandSideContribution() override;
  Real computeNeighborRightHandSideContribution() override;
  Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;
  Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;
  void setupFaceData(const FaceInfo * face_info) override;

protected:
  const RhieChowFaceFluxProvider & _mass_flux_provider;

private:
  std::pair<Real, Real> _advected_interp_coeffs;
  Real _volumetric_face_flux;

  const Moose::Functor<ADReal> * const _u_slip;
  const Moose::Functor<ADReal> * const _v_slip;
  const Moose::Functor<ADReal> * const _w_slip;

  bool _add_slip_model;
  Moose::FV::InterpMethod _advected_interp_method;
};
