//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVSRFFunctorMaterial.h"
#include "NS.h" // Variable Term Names
#include "NavierStokesMethods.h"

registerMooseObject("OpenPronghornApp", LinearFVSRFFunctorMaterial);

InputParameters
LinearFVSRFFunctorMaterial::validParams()
{
  auto params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Creates functors for the Single Rotating Frame approach");

  // Metacenter origin
  params.addRequiredParam<RealVectorValue>("mc_origin", "Vector Coordinates of the metacenter reference frame w.r.t the origin");

  // Inclination angle parameters
  params.addParam<MooseEnum>("SRF_input_mode",MooseEnum("fixed pitch_yaw_roll", "pitch_yaw_roll"),"SRF input mode");

  // Static-angle mode inputs
  params.addParam<MooseFunctorName>("pitch_angle_fixed", "Constant pitch angle in degrees");
  params.addParam<MooseFunctorName>("yaw_angle_fixed",   "Constant yaw angle in degrees");
  params.addParam<MooseFunctorName>("roll_angle_fixed",  "Constant roll angle in degrees");

  // Pitch-yaw-roll sinusoidal mode inputs
  params.addParam<MooseFunctorName>("pitch_amp", "Pitch amplitude in degrees");
  params.addParam<MooseFunctorName>("pitch_per", "Pitch period [s]");
  params.addParam<MooseFunctorName>("pitch_pha", 0.0, "Pitch phase shift in degrees");

  params.addParam<MooseFunctorName>("yaw_amp", "Yaw amplitude in degrees");
  params.addParam<MooseFunctorName>("yaw_per", "Yaw period [s]");
  params.addParam<MooseFunctorName>("yaw_pha", 0.0, "Yaw phase shift in degrees");

  params.addParam<MooseFunctorName>("roll_amp", "Roll amplitude in degrees");
  params.addParam<MooseFunctorName>("roll_per", "Roll period [s]");
  params.addParam<MooseFunctorName>("roll_pha", 0.0, "Roll phase shift in degrees");

  // Angular velocity parameters
  params.addParam<MooseFunctorName>("pitch_omega_fixed", 0.0, "Constant pitch angular velocity in rad/s");
  params.addParam<MooseFunctorName>("yaw_omega_fixed", 0.0,  "Constant yaw angular velocity in rad/s");
  params.addParam<MooseFunctorName>("roll_omega_fixed", 0.0, "Constant roll angular velocity in rad/s");

  // Angular acceleration parameters
  params.addParam<MooseFunctorName>("pitch_omegadot_fixed", 0.0, "Constant pitch angular acceleration in rad/s^2");
  params.addParam<MooseFunctorName>("yaw_omegadot_fixed", 0.0,  "Constant yaw angular acceleration in rad/s^2");
  params.addParam<MooseFunctorName>("roll_omegadot_fixed", 0.0, "Constant roll angular acceleration in rad/s^2");

  return params;
}

LinearFVSRFFunctorMaterial::LinearFVSRFFunctorMaterial(const InputParameters & params)
  : FunctorMaterial(params),
    _mc_origin(getParam<RealVectorValue>("mc_origin")),
    _SRF_input_mode(getParam<MooseEnum>("SRF_input_mode")),
    _pitch_angle_fixed(isParamValid("pitch_angle_fixed")
                            ? &getFunctor<Real>("pitch_angle_fixed")
                            : nullptr),
    _yaw_angle_fixed(isParamValid("yaw_angle_fixed")
                          ? &getFunctor<Real>("yaw_angle_fixed")
                          : nullptr),
    _roll_angle_fixed(isParamValid("roll_angle_fixed")
                           ? &getFunctor<Real>("roll_angle_fixed")
                           : nullptr),

    _pitch_omega_fixed(getFunctor<Real>("pitch_omega_fixed")),
    _yaw_omega_fixed(getFunctor<Real>("yaw_omega_fixed")),
    _roll_omega_fixed(getFunctor<Real>("roll_omega_fixed")),
    _pitch_omegadot_fixed(getFunctor<Real>("pitch_omegadot_fixed")),
    _yaw_omegadot_fixed(getFunctor<Real>("yaw_omegadot_fixed")),
    _roll_omegadot_fixed(getFunctor<Real>("roll_omegadot_fixed")),

    _pitch_amp(isParamValid("pitch_amp") ? &getFunctor<Real>("pitch_amp") : nullptr),
    _pitch_per(isParamValid("pitch_per") ? &getFunctor<Real>("pitch_per") : nullptr),
    _pitch_pha(getFunctor<Real>("pitch_pha")),
    _yaw_amp(isParamValid("yaw_amp") ? &getFunctor<Real>("yaw_amp") : nullptr),
    _yaw_per(isParamValid("yaw_per") ? &getFunctor<Real>("yaw_per") : nullptr),
    _yaw_pha(getFunctor<Real>("yaw_pha")),
    _roll_amp(isParamValid("roll_amp") ? &getFunctor<Real>("roll_amp") : nullptr),
    _roll_per(isParamValid("roll_per") ? &getFunctor<Real>("roll_per") : nullptr),
    _roll_pha(getFunctor<Real>("roll_pha"))
{
  /// Vector from Metacenter to centroid.
  addFunctorProperty<RealVectorValue>(
      "r_mc",
      [this](const auto & r, const auto & /*t*/) -> RealVectorValue
      {
        const Point xc = r.getPoint();
        return RealVectorValue(xc - _mc_origin);
      });

  if (_SRF_input_mode == "fixed")
  {
    if (!(_pitch_angle_fixed && _yaw_angle_fixed && _roll_angle_fixed))
      mooseError("For SRF_input_mode = 'fixed', provide "
                 "'pitch_angle_fixed', 'yaw_angle_fixed', and 'roll_angle_fixed'.");

    addFunctorProperty<Real>(
        "pitch_angle",
        [this](const auto & r, const auto & t) -> Real
        {
          return (*_pitch_angle_fixed)(r, t) * (libMesh::pi / 180.0);
        });

    addFunctorProperty<Real>(
        "yaw_angle",
        [this](const auto & r, const auto & t) -> Real
        {
          return (*_yaw_angle_fixed)(r, t) * (libMesh::pi / 180.0);
        });

    addFunctorProperty<Real>(
        "roll_angle",
        [this](const auto & r, const auto & t) -> Real
        {
          return (*_roll_angle_fixed)(r, t) * (libMesh::pi / 180.0);
        });

    addFunctorProperty<RealVectorValue>(
        "omega_brf",
        [this](const auto & r, const auto & t) -> RealVectorValue
        {
          return RealVectorValue(_roll_omega_fixed(r, t),_pitch_omega_fixed(r, t),_yaw_omega_fixed(r, t));
        });

    // Scalar aliases
    addFunctorProperty<Real>(
        "omega_roll",
        [this](const auto & r, const auto & t) -> Real
        {
          return _roll_omega_fixed(r, t);
        });

    addFunctorProperty<Real>(
        "omega_pitch",
        [this](const auto & r, const auto & t) -> Real
        {
          return _pitch_omega_fixed(r, t);
        });

    addFunctorProperty<Real>(
        "omega_yaw",
        [this](const auto & r, const auto & t) -> Real
        {
          return _yaw_omega_fixed(r, t);
        });

    addFunctorProperty<RealVectorValue>(
        "omega_dot_brf",
        [this](const auto & r, const auto & t) -> RealVectorValue
        {
          return RealVectorValue(_roll_omegadot_fixed(r, t),_pitch_omegadot_fixed(r, t),_yaw_omegadot_fixed(r, t));
        });

    addFunctorProperty<Real>(
        "omega_dot_roll",
        [this](const auto & r, const auto & t) -> Real
        {
          return _roll_omegadot_fixed(r, t);
        });

    addFunctorProperty<Real>(
        "omega_dot_pitch",
        [this](const auto & r, const auto & t) -> Real
        {
          return _pitch_omegadot_fixed(r, t);
        });

    addFunctorProperty<Real>(
        "omega_dot_yaw",
        [this](const auto & r, const auto & t) -> Real
        {
          return _yaw_omegadot_fixed(r, t);
        });
  }
  else if (_SRF_input_mode == "pitch_yaw_roll")
  {
    if (!(_pitch_amp && _pitch_per && _yaw_amp && _yaw_per && _roll_amp && _roll_per))
      mooseError("For SRF_input_mode = 'pitch_yaw_roll', provide "
                 "'pitch_amp', 'pitch_per', 'yaw_amp', 'yaw_per', 'roll_amp', and 'roll_per'. "
                 "The phase parameters are optional.");

    const Real deg_to_rad = libMesh::pi / 180.0;

    auto eval_angle =
        [this, deg_to_rad](const Moose::Functor<Real> & amp_deg,
                           const Moose::Functor<Real> & per,
                           const Moose::Functor<Real> & pha_deg,
                           const std::string & label,
                           const auto & r,
                           const auto & t) -> Real
    {
      using std::sin;

      const Real T = per(r, t);
      if (T <= 0.0)
        mooseError(label, "_per must be positive.");

      const Real A_rad = amp_deg(r, t) * deg_to_rad;
      const Real phase_rad = pha_deg(r, t) * deg_to_rad;
      const Real w = 2.0 * libMesh::pi / T;

      return A_rad * sin(w * _t + phase_rad);
    };

    auto eval_rate =
        [this, deg_to_rad](const Moose::Functor<Real> & amp_deg,
                           const Moose::Functor<Real> & per,
                           const Moose::Functor<Real> & pha_deg,
                           const std::string & label,
                           const auto & r,
                           const auto & t) -> Real
    {
      using std::cos;

      const Real T = per(r, t);
      if (T <= 0.0)
        mooseError(label, "_per must be positive.");

      const Real A_rad = amp_deg(r, t) * deg_to_rad;
      const Real phase_rad = pha_deg(r, t) * deg_to_rad;
      const Real w = 2.0 * libMesh::pi / T;

      return A_rad * w * cos(w * _t + phase_rad);
    };

    auto eval_accel =
        [this, deg_to_rad](const Moose::Functor<Real> & amp_deg,
                           const Moose::Functor<Real> & per,
                           const Moose::Functor<Real> & pha_deg,
                           const std::string & label,
                           const auto & r,
                           const auto & t) -> Real
    {
      using std::sin;

      const Real T = per(r, t);
      if (T <= 0.0)
        mooseError(label, "_per must be positive.");

      const Real A_rad = amp_deg(r, t) * deg_to_rad;
      const Real phase_rad = pha_deg(r, t) * deg_to_rad;
      const Real w = 2.0 * libMesh::pi / T;

      return -A_rad * w * w * sin(w * _t + phase_rad);
    };

    // Angles
    addFunctorProperty<Real>(
        "pitch_angle",
        [this, eval_angle](const auto & r, const auto & t) -> Real
        {
          return eval_angle(*_pitch_amp, *_pitch_per, _pitch_pha, "pitch", r, t);
        });

    addFunctorProperty<Real>(
        "yaw_angle",
        [this, eval_angle](const auto & r, const auto & t) -> Real
        {
          return eval_angle(*_yaw_amp, *_yaw_per, _yaw_pha, "yaw", r, t);
        });

    addFunctorProperty<Real>(
        "roll_angle",
        [this, eval_angle](const auto & r, const auto & t) -> Real
        {
          return eval_angle(*_roll_amp, *_roll_per, _roll_pha, "roll", r, t);
        });

    // Body-frame omega vector
    auto omega_brf_value =
        [this, eval_angle, eval_rate](const auto & r, const auto & t) -> RealVectorValue
    {
      using std::cos;
      using std::sin;

      const Real pitch = eval_angle(*_pitch_amp, *_pitch_per, _pitch_pha, "pitch", r, t);
      const Real yaw   = eval_angle(*_yaw_amp, *_yaw_per, _yaw_pha, "yaw", r, t);
      const Real roll  = eval_angle(*_roll_amp, *_roll_per, _roll_pha, "roll", r, t);

      const Real pitch_dot = eval_rate(*_pitch_amp, *_pitch_per, _pitch_pha, "pitch", r, t);
      const Real yaw_dot   = eval_rate(*_yaw_amp, *_yaw_per, _yaw_pha, "yaw", r, t);
      const Real roll_dot  = eval_rate(*_roll_amp, *_roll_per, _roll_pha, "roll", r, t);

      // Body-frame angular velocity vector
      const RealVectorValue omega_brf(
          roll_dot + sin(yaw) * pitch_dot,
          cos(roll) * cos(yaw) * pitch_dot + sin(roll) * yaw_dot,
          -sin(roll) * cos(yaw) * pitch_dot + cos(roll) * yaw_dot);

      return omega_brf;
    };

    // Body-frame angular acceleration vector
    auto omegadot_brf_value =
        [this, eval_angle, eval_rate, eval_accel](const auto & r, const auto & t) -> RealVectorValue
    {
      using std::cos;
      using std::sin;

      const Real pitch = eval_angle(*_pitch_amp, *_pitch_per, _pitch_pha, "pitch", r, t);
      const Real yaw   = eval_angle(*_yaw_amp, *_yaw_per, _yaw_pha, "yaw", r, t);
      const Real roll  = eval_angle(*_roll_amp, *_roll_per, _roll_pha, "roll", r, t);

      const Real pitch_dot = eval_rate(*_pitch_amp, *_pitch_per, _pitch_pha, "pitch", r, t);
      const Real yaw_dot   = eval_rate(*_yaw_amp, *_yaw_per, _yaw_pha, "yaw", r, t);
      const Real roll_dot  = eval_rate(*_roll_amp, *_roll_per, _roll_pha, "roll", r, t);

      const Real pitch_ddot = eval_accel(*_pitch_amp, *_pitch_per, _pitch_pha, "pitch", r, t);
      const Real yaw_ddot   = eval_accel(*_yaw_amp, *_yaw_per, _yaw_pha, "yaw", r, t);
      const Real roll_ddot  = eval_accel(*_roll_amp, *_roll_per, _roll_pha, "roll", r, t);

      // Body-frame angular acceleration vector
      const RealVectorValue omegadot_brf(
          roll_ddot + cos(yaw) * yaw_dot * pitch_dot + sin(yaw) * pitch_ddot,
          cos(roll) * cos(yaw) * pitch_ddot + sin(roll) * yaw_ddot +
              cos(roll) * roll_dot * yaw_dot - sin(roll) * cos(yaw) * roll_dot * pitch_dot -
              cos(roll) * sin(yaw) * yaw_dot * pitch_dot,
          -sin(roll) * cos(yaw) * pitch_ddot + cos(roll) * yaw_ddot -
              sin(roll) * roll_dot * yaw_dot - cos(roll) * cos(yaw) * roll_dot * pitch_dot +
              sin(roll) * sin(yaw) * yaw_dot * pitch_dot);

      return omegadot_brf;
    };

    addFunctorProperty<RealVectorValue>(
        "omega_brf",
        [this, omega_brf_value](const auto & r, const auto & t) -> RealVectorValue
        {
          return omega_brf_value(r, t);
        });

    // Scalar aliases from body-frame vectors: (x,y,z) = (roll,pitch,yaw)
    addFunctorProperty<Real>(
        "omega_roll",
        [this, omega_brf_value](const auto & r, const auto & t) -> Real
        {
          return omega_brf_value(r, t)(0);
        });

    addFunctorProperty<Real>(
        "omega_pitch",
        [this, omega_brf_value](const auto & r, const auto & t) -> Real
        {
          return omega_brf_value(r, t)(1);
        });

    addFunctorProperty<Real>(
        "omega_yaw",
        [this, omega_brf_value](const auto & r, const auto & t) -> Real
        {
          return omega_brf_value(r, t)(2);
        });

    addFunctorProperty<RealVectorValue>(
        "omega_dot_brf",
        [this, omegadot_brf_value](const auto & r, const auto & t) -> RealVectorValue
        {
          return omegadot_brf_value(r, t);
        });

    addFunctorProperty<Real>(
        "omega_dot_roll",
        [this, omegadot_brf_value](const auto & r, const auto & t) -> Real
        {
          return omegadot_brf_value(r, t)(0);
        });

    addFunctorProperty<Real>(
        "omega_dot_pitch",
        [this, omegadot_brf_value](const auto & r, const auto & t) -> Real
        {
          return omegadot_brf_value(r, t)(1);
        });

    addFunctorProperty<Real>(
        "omega_dot_yaw",
        [this, omegadot_brf_value](const auto & r, const auto & t) -> Real
        {
          return omegadot_brf_value(r, t)(2);
        });
  }
}
