//* This file is part of OpenPronghorn.
//* https://github.com/idaholab/open_pronghorn
//* https://mooseframework.inl.gov/open_pronghorn
//*
//* OpenPronghorn is powered by the MOOSE Framework
//* https://mooseframework.inl.gov
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Copyright 2025, Battelle Energy Alliance, LLC
//* ALL RIGHTS RESERVED
//*

#pragma once

#include "WCNSLinearFVTurbulencePhysics.h"

/**
 * Creates the objects needed for a k-epsilon turbulence model with a linear finite volume
 * discretization, using the OpenPronghorn k-epsilon objects. With regard to the MOOSE
 * WCNSLinearFVTurbulencePhysics, this adds the two-layer, low-Reynolds and realizable k-epsilon
 * variants, along with the buoyancy, compressibility, Yap, non-linear constitutive and
 * curvature corrections.
 */
class AdvancedKEpsilonPhysics final : public WCNSLinearFVTurbulencePhysics
{
public:
  static InputParameters validParams();

  AdvancedKEpsilonPhysics(const InputParameters & parameters);

protected:
  virtual void addAuxiliaryVariables() override;
  virtual void addAuxiliaryKernels() override;
  virtual void addTurbulentViscosityAuxKernel() override;
  virtual void addKEpsilonSink() override;

private:
  /// Whether the selected variant and corrections require a distance to the nearest wall
  bool needsWallDistance() const;

  /**
   * The OpenPronghorn k-epsilon objects take several closure coefficients as constants, while
   * the Physics base class declares them as functors. Converts such a parameter, erroring out
   * if the user provided a non-constant functor.
   * @param param_name name of the functor parameter to convert
   */
  Real getConstantCoefficient(const std::string & param_name) const;

  /// Parameters forwarded as-is to the OpenPronghorn k-epsilon objects. Only the parameters
  /// the user actually set are forwarded, so each object keeps its own defaults and can
  /// detect the model options that do not apply to it
  static const std::vector<std::string> _variant_parameters;

  /// Name of the functor providing the distance to the nearest wall
  const MooseFunctorName _wall_distance_name;

  /// Whether this Physics creates the wall distance auxiliary variable and kernel
  const bool _create_wall_distance;
};
