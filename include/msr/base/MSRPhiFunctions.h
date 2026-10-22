#pragma once

#include "MooseTypes.h"

#include <utility>
#include <vector>

namespace MSR
{

using DenseMatrix = std::vector<std::vector<Real>>;

/**
 * Dense phi-function actions for EPI3V.
 *
 * Small systems use a fixed contiguous scaling-and-squaring implementation.
 * Larger systems use the general dense implementation.
 */
class MSRPhiFunctions
{
public:
  static std::vector<Real> phi1Action(const DenseMatrix & A,
                                      const std::vector<Real> & v);

  static std::vector<Real> phi3Action(const DenseMatrix & A,
                                      const std::vector<Real> & v);

  static std::pair<std::vector<Real>, std::vector<Real>>
  phi1Phi3Actions(const DenseMatrix & A,
                  const std::vector<Real> & v1,
                  const std::vector<Real> & v3);

private:
  static std::vector<Real> phiAction(const DenseMatrix & A,
                                     const std::vector<Real> & v,
                                     unsigned int k);

  static std::vector<Real> matVec(const DenseMatrix & A,
                                  const std::vector<Real> & v);
};

} // namespace MSR
