//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WallDistanceAux.h"

#include "libmesh/parallel_algebra.h"

registerMooseObject("OpenPronghornApp", WallDistanceAux);

InputParameters
WallDistanceAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Computes the distance from the nearest wall.");
  params.addRequiredParam<std::vector<BoundaryName>>("walls",
                                                     "Boundaries that correspond to solid walls.");
  return params;
}

WallDistanceAux::WallDistanceAux(const InputParameters & parameters)
  : AuxKernel(parameters), _wall_boundary_names(getParam<std::vector<BoundaryName>>("walls"))
{
  if (!dynamic_cast<MooseVariableFV<Real> *>(&_var))
    paramError("variable",
               "'",
               name(),
               "' is currently programmed to use finite volume machinery, so make sure that '",
               _var.name(),
               "' is a finite volume variable.");
}

void
WallDistanceAux::initialSetup()
{
  meshChanged();
}

void
WallDistanceAux::meshChanged()
{
  // Get the ids of the wall boundaries
  std::vector<BoundaryID> vec_ids = _mesh.getBoundaryIDs(_wall_boundary_names, true);
  std::vector<BoundaryName> invalid_boundaries;
  for (const auto i : index_range(vec_ids))
    if (vec_ids[i] == BoundaryInfo::invalid_id)
      invalid_boundaries.push_back(_wall_boundary_names[i]);
  if (!invalid_boundaries.empty())
    paramError("walls",
               "The following boundaries do not exist in the mesh: ",
               Moose::stringify(invalid_boundaries));

  // Map of boundary IDs to element ids
  const auto & bnd_to_elem_map = _mesh.getBoundariesToActiveSemiLocalElemIds();

  // Use a local set during collection to deduplicate faces visited from both elem and neighbor
  // sides
  std::unordered_set<Point> local_set;

  // Loop through wall boundaries
  for (const auto & bid : vec_ids)
  {
    // Get element IDs associated with boundary
    auto search = bnd_to_elem_map.find(bid);
    // If boundary doesn't exist locally, just continue
    if (search == bnd_to_elem_map.end())
      continue;
    const auto & bnd_elem_ids = search->second;

    // Loop through boundary elements to gather face centroids
    for (const auto bnd_elem_id : bnd_elem_ids)
    {
      // Get a pointer to the element
      const auto elem = _mesh.elemPtr(bnd_elem_id);
      // This shouldn't happen
      mooseAssert(elem,
                  "Boundary element " + std::to_string(bnd_elem_id) + " does not exist on mesh.");

      // Get side(s) where boundary is
      const auto sides = _mesh.getMesh().get_boundary_info().sides_with_boundary_id(elem, bid);
      // This shouldn't happen
      mooseAssert(sides.size() > 0,
                  "Element " + std::to_string(bnd_elem_id) + " does not contain boundary " +
                      std::to_string(bid));

      // Loop through sides and collect face centroid
      for (const auto side : sides)
      {
        auto fi = _mesh.faceInfo(elem, side);
        // It's possible that we are on an internal boundary
        if (!fi)
        {
          const Elem * const neigh = elem->neighbor_ptr(side);
          mooseAssert(neigh,
                      "In WallDistanceAux, we could not find a face information object with elem "
                      "and side, and we are on an external boundary. This shouldn't happen.");
          const auto neigh_side = neigh->which_neighbor_am_i(elem);
          fi = _mesh.faceInfo(neigh, neigh_side);
          mooseAssert(fi, "We should have a face info for either the elem or neigh side");
        }

        // Collect centroid
        local_set.insert(fi->faceCentroid());
      }
    }
  }

  // If the mesh is distributed, we have to allgather
  if (_mesh.isDistributedMesh())
    _communicator.set_union(local_set);

  if (local_set.empty())
    paramError("walls", "Could not find boundaries");

  _boundary_points.assign(local_set.begin(), local_set.end());
  _kd_tree = std::make_unique<KDTree>(_boundary_points, /*max_leaf_size=*/10);
}

Real
WallDistanceAux::computeValue()
{
  mooseAssert(_kd_tree, "WallDistanceAux has not been initialized yet.");

  std::vector<std::size_t> idx(1);
  std::vector<Real> d2(1);
  _kd_tree->neighborSearch(_q_point[_qp], 1, idx, d2);
  return std::sqrt(d2[0]);
}
