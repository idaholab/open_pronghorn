//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MathFVUtils.h"
#include "LinearFVFluxKernel.h"

class LinearFVBoundaryCondition;

/**
 * Linear finite-volume flux kernel that adds the aerosol drift-flux contribution
 * to a transported scalar (e.g. aerosol mass or volume fraction).
 *
 * The conservative equation for rho * C is
 *
 *   d(rho C)/dt + div( rho U C ) + div( rho U_drift C ) = ...
 *
 * This kernel provides the flux term:
 *
 *   J_drift . n = (rho_face * U_drift_face . n) * C_face
 *
 * where U_drift is the drift velocity resulting from gravitational settling and
 * thermophoresis, computed internally from local thermofluid and particle properties.
 */
class LinearFVAerosolDriftFlux : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();

  LinearFVAerosolDriftFlux(const InputParameters & params);

  /// Matrix contribution on the element side of an internal face
  virtual Real computeElemMatrixContribution() override;

  /// Matrix contribution on the neighbor side of an internal face
  virtual Real computeNeighborMatrixContribution() override;

  /// RHS contribution on the element side of an internal face
  virtual Real computeElemRightHandSideContribution() override;

  /// RHS contribution on the neighbor side of an internal face
  virtual Real computeNeighborRightHandSideContribution() override;

  /// For now we do not add any explicit boundary matrix terms here
  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition &) override
  {
    return 0.0;
  }

  /// Boundary RHS contribution (currently zero, drift handled via matrix flux)
  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;

  /**
   * Set the current FaceInfo object and cache drift mass flux + interpolation
   * coefficients. We override this so that face work is done only once per face.
   */
  virtual void setupFaceData(const FaceInfo * face_info) override;

protected:
  /// Compute drift mass flux on the current face: rho_face * (U_drift . n)
  void computeFlux();

  /// Compute gravitational settling drift velocity at the face
  RealVectorValue
  computeGravitationalDrift(const Moose::FaceArg & face_arg, Real mu_face, Real rho_face) const;

  /// Compute thermophoretic drift velocity at the face (if enabled)
  RealVectorValue
  computeThermophoreticDrift(const Moose::FaceArg & face_arg, Real mu_face, Real rho_face) const;

  /// Spatial dimension of the simulation
  const unsigned int _dim;

  /// Carrier fluid density (for rho * C conservative form)
  const Moose::Functor<Real> & _rho;

  /// Carrier fluid dynamic viscosity
  const Moose::Functor<Real> & _mu;

  /// Carrier fluid temperature
  const Moose::Functor<Real> & _T;

  /// Interpolation method for density at faces
  const Moose::FV::InterpMethod _density_interp_method;

  /// Interpolation method for viscosity at faces
  const Moose::FV::InterpMethod _viscosity_interp_method;

  /// Particle diameter [m]
  const Real _d_p;

  /// Particle density [kg/m^3]
  const Real _rho_p;

  /// Mean free path of the carrier gas [m]
  const Real _lambda_mfp;

  /// Optional user-specified thermophoretic coefficient (dimensionless)
  const Real _k_th_const;

  /// If true, use the user-specified constant k_th; otherwise use a simple
  /// diameter-based correlation (small-particle approximation)
  const bool _use_constant_kth;

  /// Gravitational acceleration vector [m/s^2]
  const RealVectorValue _gravity;

  /// Face drift mass flux: m_drift = rho_face * (U_drift . n)
  Real _face_mass_flux;

  /// Upwind/central interpolation coefficients for the scalar at the face
  std::pair<Real, Real> _interp_coeffs;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;
};
