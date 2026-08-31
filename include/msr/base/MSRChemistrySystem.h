#pragma once

#include "MoltenSaltRadiolysisData.h"

#include <map>
#include <string>
#include <vector>

namespace MSR
{
class MSRChemistrySystem
{
public:
  MSRChemistrySystem(const std::vector<std::string> & species,
                     const std::vector<ReactionData> & reactions);

  MSRChemistrySystem(const std::vector<std::string> & species,
                     const std::vector<ReactionData> & reactions,
                     const std::vector<Real> & constant_sources);

  std::vector<Real> evaluateRHS(const std::vector<Real> & y, Real temperature) const;

  std::vector<std::vector<Real>>
  evaluateJacobian(const std::vector<Real> & y, Real temperature) const;

  const std::vector<std::string> & species() const { return _species; }
  const std::vector<Real> & constantSources() const { return _constant_sources; }

private:
  std::vector<std::string> _species;
  std::vector<ReactionData> _reactions;
  std::vector<Real> _constant_sources;
  std::map<std::string, std::size_t> _species_index;
};
} // namespace MSR
