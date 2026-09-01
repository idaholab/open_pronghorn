#pragma once

#include "TimeIntegrator.h"
#include "MSRChemistrySystem.h"
#include "MoltenSaltRadiolysisData.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Cell-local EPI3V time integrator for molten-salt radiolysis chemistry.
 *
 * Temperature and dose rate may be supplied as constants or cell-local FV variables.
 */
class MSREPI3VTimeIntegrator : public TimeIntegrator
{
public:
  static InputParameters validParams();

  MSREPI3VTimeIntegrator(const InputParameters & parameters);

  virtual int order() override { return 3; }

  virtual bool isExplicit() const override { return true; }

  virtual bool overridesSolve() const override { return true; }

  virtual void solve() override;

  virtual void computeTimeDerivatives() override;

  void computeADTimeDerivatives(ADReal & ad_u_dot,
                                const dof_id_type & dof,
                                ADReal & ad_u_dotdot) const override;

  virtual void postResidual(NumericVector<Number> & residual) override;

  Real lastErrorNorm() const { return _last_error_norm; }

protected:
  std::vector<std::string> _species;
  std::vector<MSR::ReactionData> _reactions;
  std::unique_ptr<MSR::MSRChemistrySystem> _chemistry;

  /// Constant chemistry temperature [K].
  const Real _temperature;

  /// Optional cell-local FV temperature variable. When nonempty, overrides _temperature.
  const std::string _temperature_variable;

  /// Constant volumetric dose rate [J/m^3/s].
  const Real _dose_rate;

  /// Optional cell-local FV dose-rate variable. When nonempty, overrides _dose_rate.
  const std::string _dose_rate_variable;

  /// Radiolytic source coefficient per unit dose rate for each species.
  std::vector<Real> _source_per_unit_dose;

  const Real _rtol;
  const Real _atol;

  Real _last_error_norm;
};
