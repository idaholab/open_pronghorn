//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"
#include "libmesh/vector_value.h"

namespace NS
{
  namespace SRF
  {

    using libMesh::RealVectorValue;

    RealVectorValue
    rotateVectorInertialToBody(const RealVectorValue & v_inertial,
                              const Real & pitch,
                              const Real & yaw,
                              const Real & roll);

    RealVectorValue
    rotateVectorBodyToInertial(const RealVectorValue & v_body,
                              const Real & pitch,
                              const Real & yaw,
                              const Real & roll);

  }
}
