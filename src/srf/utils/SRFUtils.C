//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SRFUtils.h"

#include <cmath>

namespace NS
{
  namespace SRF
  {

    RealVectorValue
    rotateVectorInertialToBody(const RealVectorValue & v_inertial,
                              const Real & pitch,
                              const Real & yaw,
                              const Real & roll)
    {
      using std::cos;
      using std::sin;

      const Real cp = cos(pitch);
      const Real sp = sin(pitch);
      const Real cy = cos(yaw);
      const Real sy = sin(yaw);
      const Real cr = cos(roll);
      const Real sr = sin(roll);

      return RealVectorValue((cy * cp) * v_inertial(0) + sy * v_inertial(1) - cy * sp * v_inertial(2),
                            (sr * sp - cr * sy * cp) * v_inertial(0) + cr * cy * v_inertial(1) +
                                (sr * cp + cr * sy * sp) * v_inertial(2),
                            (cr * sp + sr * sy * cp) * v_inertial(0) - sr * cy * v_inertial(1) +
                                (cr * cp - sr * sy * sp) * v_inertial(2));
    }

    RealVectorValue
    rotateVectorBodyToInertial(const RealVectorValue & v_body,
                              const Real & pitch,
                              const Real & yaw,
                              const Real & roll)
    {
      using std::cos;
      using std::sin;

      const Real cp = cos(pitch);
      const Real sp = sin(pitch);
      const Real cy = cos(yaw);
      const Real sy = sin(yaw);
      const Real cr = cos(roll);
      const Real sr = sin(roll);

      return RealVectorValue(
          (cy * cp) * v_body(0) + (sr * sp - cr * sy * cp) * v_body(1) +
              (cr * sp + sr * sy * cp) * v_body(2),

          sy * v_body(0) + (cr * cy) * v_body(1) - (sr * cy) * v_body(2),

          (-cy * sp) * v_body(0) + (sr * cp + cr * sy * sp) * v_body(1) +
              (cr * cp - sr * sy * sp) * v_body(2));
    }

  }
}
