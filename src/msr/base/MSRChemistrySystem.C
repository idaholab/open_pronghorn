#include "MSRChemistrySystem.h"
#include "MooseError.h"

#include <cmath>

namespace MSR
{
MSRChemistrySystem::MSRChemistrySystem(const std::vector<std::string> & species,
                                       const std::vector<ReactionData> & reactions)
  : MSRChemistrySystem(species, reactions, std::vector<Real>(species.size(), 0.0))
{
}

MSRChemistrySystem::MSRChemistrySystem(const std::vector<std::string> & species,
                                       const std::vector<ReactionData> & reactions,
                                       const std::vector<Real> & constant_sources)
  : _species(species), _reactions(reactions), _constant_sources(constant_sources)
{
  if (_constant_sources.size() != _species.size())
    mooseError("MSRChemistrySystem: constant source vector has size ",
               _constant_sources.size(),
               " but the chemistry system contains ",
               _species.size(),
               " species.");

  for (std::size_t i = 0; i < _species.size(); ++i)
  {
    const auto inserted = _species_index.emplace(_species[i], i);
    if (!inserted.second)
      mooseError("MSRChemistrySystem: duplicate species '", _species[i], "'.");
  }
}

std::vector<Real>
MSRChemistrySystem::evaluateRHS(const std::vector<Real> & y, Real temperature) const
{
  if (y.size() != _species.size())
    mooseError("MSRChemistrySystem: concentration vector has size ",
               y.size(),
               " but the chemistry system contains ",
               _species.size(),
               " species.");

  std::vector<Real> rhs = _constant_sources;

  for (const auto & reaction : _reactions)
  {
    Real reaction_rate = reaction.rate(temperature);

    for (const auto & reactant : reaction.reactants)
    {
      const auto it = _species_index.find(reactant.species);
      if (it == _species_index.end())
      {
        reaction_rate = 0.0;
        break;
      }

      const Real concentration = y[it->second];
      if (concentration < 0.0)
        mooseError("MSRChemistrySystem: negative concentration for species '",
                   reactant.species,
                   "' while evaluating reaction '",
                   reaction.name,
                   "'.");

      reaction_rate *= std::pow(concentration, reactant.order);
    }

    for (const auto & reactant : reaction.reactants)
    {
      const auto it = _species_index.find(reactant.species);
      if (it != _species_index.end())
        rhs[it->second] -= reactant.coeff * reaction_rate;
    }

    for (const auto & product : reaction.products)
    {
      const auto it = _species_index.find(product.species);
      if (it != _species_index.end())
        rhs[it->second] += product.coeff * reaction_rate;
    }
  }

  return rhs;
}

std::vector<std::vector<Real>>
MSRChemistrySystem::evaluateJacobian(const std::vector<Real> & y, Real temperature) const
{
  if (y.size() != _species.size())
    mooseError("MSRChemistrySystem: concentration vector has size ",
               y.size(),
               " but the chemistry system contains ",
               _species.size(),
               " species.");

  const std::size_t n = _species.size();
  std::vector<std::vector<Real>> jacobian(n, std::vector<Real>(n, 0.0));

  for (const auto & reaction : _reactions)
  {
    const Real k = reaction.rate(temperature);

    for (const auto & differentiated_reactant : reaction.reactants)
    {
      const auto differentiated_it = _species_index.find(differentiated_reactant.species);
      if (differentiated_it == _species_index.end())
        continue;

      const std::size_t column = differentiated_it->second;
      Real dr_dy = k * differentiated_reactant.order;

      for (const auto & reactant : reaction.reactants)
      {
        const auto it = _species_index.find(reactant.species);
        if (it == _species_index.end())
        {
          dr_dy = 0.0;
          break;
        }

        const Real concentration = y[it->second];
        if (concentration < 0.0)
          mooseError("MSRChemistrySystem: negative concentration for species '",
                     reactant.species,
                     "' while evaluating reaction '",
                     reaction.name,
                     "'.");

        Real exponent = reactant.order;
        if (reactant.species == differentiated_reactant.species)
          exponent -= 1.0;

        if (exponent == 0.0)
          continue;

        dr_dy *= std::pow(concentration, exponent);
      }

      for (const auto & reactant : reaction.reactants)
      {
        const auto it = _species_index.find(reactant.species);
        if (it != _species_index.end())
          jacobian[it->second][column] -= reactant.coeff * dr_dy;
      }

      for (const auto & product : reaction.products)
      {
        const auto it = _species_index.find(product.species);
        if (it != _species_index.end())
          jacobian[it->second][column] += product.coeff * dr_dy;
      }
    }
  }

  return jacobian;
}
} // namespace MSR
