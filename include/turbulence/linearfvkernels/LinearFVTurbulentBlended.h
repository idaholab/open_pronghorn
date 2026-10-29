//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVTurbulentAdvection.h"
#include "LinearFVTurbulentDiffusion.h"
#include "LinearFVTurbulentTimeDerivative.h"
#include "LinearFVFluxKernel.h"
#include "LinearFVElementalKernel.h"
#include "NavierStokesMethods.h"
#include "TurbulenceMethods.h"

#include <optional>
#include <type_traits>
#include <unordered_set>

/**
 * Non-template interface shared by the blended turbulence kernels, so other objects (e.g. the
 * two-layer epsilon source/sink) can identify them and check which blending factor they use.
 */
class LinearFVTurbulentBlendedInterface
{
public:
  LinearFVTurbulentBlendedInterface(const InputParameters & params)
    : _blending_name(
          params.isParamValid("blending_factor")
              ? std::optional<MooseFunctorName>(params.get<MooseFunctorName>("blending_factor"))
              : std::nullopt)
  {
  }

  virtual ~LinearFVTurbulentBlendedInterface() = default;

  /// Name of the blending factor functor, empty if none was provided
  const std::optional<MooseFunctorName> & blendingFactorName() const { return _blending_name; }

private:
  const std::optional<MooseFunctorName> _blending_name;
};

/**
 * Turbulence kernel (advection, diffusion or time derivative) whose element rows are scaled by a
 * continuous per-element weight instead of the binary wall-bounded skip of the base kernel.
 * The face (or element) math is entirely taken from the base class.
 */
template <class Base>
class LinearFVTurbulentBlended : public Base, public LinearFVTurbulentBlendedInterface
{
public:
  static InputParameters validParams();

  LinearFVTurbulentBlended(const InputParameters & params);

  virtual void initialSetup() override;

  virtual void addMatrixContribution() override;

  virtual void addRightHandSideContribution() override;

protected:
  /// Weight of the rows of an element
  Real weight(const Elem * elem) const;

  /// Face loops of the flux kernels, with each side scaled by the weight of its element
  void addFluxMatrixContribution();
  void addFluxRightHandSideContribution();

  /// Optional blending factor
  const Moose::Functor<Real> * _blending;

  /// Wall boundaries
  const std::vector<BoundaryName> & _blend_walls;

  /// Wall-bounded elements
  std::unordered_set<const Elem *> _blend_wall_bounded;
};

template <class Base>
InputParameters
LinearFVTurbulentBlended<Base>::validParams()
{
  InputParameters params = Base::validParams();
  params.addClassDescription(
      "Turbulence kernel whose rows are scaled by a blending factor and removed in wall-bounded "
      "cells.");
  params.addParam<MooseFunctorName>(
      "blending_factor",
      "Functor in [0, 1] that scales the rows of this kernel (0 removes the contribution of the "
      "kernel, 1 keeps it). Wall-bounded cells are always removed.");
  return params;
}

template <class Base>
LinearFVTurbulentBlended<Base>::LinearFVTurbulentBlended(const InputParameters & params)
  : Base(params),
    LinearFVTurbulentBlendedInterface(params),
    _blending(params.isParamValid("blending_factor")
                  ? &this->template getFunctor<Real>("blending_factor")
                  : nullptr),
    _blend_walls(this->template getParam<std::vector<BoundaryName>>("walls"))
{
}

template <class Base>
void
LinearFVTurbulentBlended<Base>::initialSetup()
{
  Base::initialSetup();
  NS::getWallBoundedElements(
      _blend_walls, this->_fe_problem, this->_subproblem, this->blockIDs(), _blend_wall_bounded);
}

template <class Base>
Real
LinearFVTurbulentBlended<Base>::weight(const Elem * elem) const
{
  return NS::turbulentBlendWeight(
      this->makeElemArg(elem), _blend_wall_bounded, _blending, this->determineState());
}

template <class Base>
void
LinearFVTurbulentBlended<Base>::addMatrixContribution()
{
  if constexpr (std::is_base_of_v<LinearFVElementalKernel, Base>)
  {
    // Elemental kernels only contribute to the diagonal of the current element
    const Real w = weight(this->_current_elem_info->elem());
    if (w == 0.0)
      return;
    const auto mx_contrib = w * this->computeMatrixContribution();
    for (auto & matrix : this->_matrices)
      (*matrix).add(this->_dof_id, this->_dof_id, mx_contrib);
  }
  else
    addFluxMatrixContribution();
}

template <class Base>
void
LinearFVTurbulentBlended<Base>::addRightHandSideContribution()
{
  if constexpr (std::is_base_of_v<LinearFVElementalKernel, Base>)
  {
    const Real w = weight(this->_current_elem_info->elem());
    if (w == 0.0)
      return;
    const auto rhs_contrib = w * this->computeRightHandSideContribution();
    for (auto & vector : this->_vectors)
      (*vector).add(this->_dof_id, rhs_contrib);
  }
  else
    addFluxRightHandSideContribution();
}

template <class Base>
void
LinearFVTurbulentBlended<Base>::addFluxMatrixContribution()
{
  const auto & fi = *this->_current_face_info;

  if (this->_current_face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    this->_dof_indices(0) = fi.elemInfo()->dofIndices()[this->_sys_num][this->_var_num];
    this->_dof_indices(1) = fi.neighborInfo()->dofIndices()[this->_sys_num][this->_var_num];

    const auto elem_matrix_contribution = this->computeElemMatrixContribution();
    const auto neighbor_matrix_contribution = this->computeNeighborMatrixContribution();

    if (this->hasBlocks(fi.elemInfo()->subdomain_id()))
    {
      const Real w = weight(fi.elemPtr());
      this->_matrix_contribution(0, 0) = w * elem_matrix_contribution;
      this->_matrix_contribution(0, 1) = w * neighbor_matrix_contribution;
    }

    if (this->hasBlocks(fi.neighborInfo()->subdomain_id()))
    {
      const Real w = weight(fi.neighborPtr());
      this->_matrix_contribution(1, 0) = -w * elem_matrix_contribution;
      this->_matrix_contribution(1, 1) = -w * neighbor_matrix_contribution;
    }
    (*this->_linear_system.matrix)
        .add_matrix(this->_matrix_contribution, this->_dof_indices.get_values());
  }
  else if (this->_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
           this->_current_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR)
  {
    mooseAssert(fi.boundaryIDs().size() == 1, "We should only have one boundary on every face.");

    LinearFVBoundaryCondition * bc_pointer =
        this->_var.getBoundaryCondition(*fi.boundaryIDs().begin());

    if (bc_pointer || this->_force_boundary_execution)
    {
      if (bc_pointer)
        bc_pointer->setupFaceData(&fi, this->_current_face_type);
      const auto matrix_contribution = this->computeBoundaryMatrixContribution(*bc_pointer);

      // Internal (for the mesh) boundaries are allowed, so the owning side matters
      const bool elem_side = this->_current_face_type == FaceInfo::VarFaceNeighbors::ELEM;
      const auto & info = elem_side ? *fi.elemInfo() : *fi.neighborInfo();
      const Real w = weight(elem_side ? fi.elemPtr() : fi.neighborPtr());
      const auto dof = info.dofIndices()[this->_sys_num][this->_var_num];
      (*this->_linear_system.matrix).add(dof, dof, w * matrix_contribution);
    }
  }
}

template <class Base>
void
LinearFVTurbulentBlended<Base>::addFluxRightHandSideContribution()
{
  const auto & fi = *this->_current_face_info;

  if (this->_current_face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    this->_dof_indices(0) = fi.elemInfo()->dofIndices()[this->_sys_num][this->_var_num];
    this->_dof_indices(1) = fi.neighborInfo()->dofIndices()[this->_sys_num][this->_var_num];

    const auto elem_rhs_contribution = this->computeElemRightHandSideContribution();
    const auto neighbor_rhs_contribution = this->computeNeighborRightHandSideContribution();

    if (this->hasBlocks(fi.elemInfo()->subdomain_id()))
      this->_rhs_contribution(0) = weight(fi.elemPtr()) * elem_rhs_contribution;
    if (this->hasBlocks(fi.neighborInfo()->subdomain_id()))
      this->_rhs_contribution(1) = weight(fi.neighborPtr()) * neighbor_rhs_contribution;

    (*this->_linear_system.rhs)
        .add_vector(this->_rhs_contribution.get_values().data(), this->_dof_indices.get_values());
  }
  else if (this->_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
           this->_current_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR)
  {
    mooseAssert(fi.boundaryIDs().size() == 1, "We should only have one boundary on every face.");

    LinearFVBoundaryCondition * bc_pointer =
        this->_var.getBoundaryCondition(*fi.boundaryIDs().begin());

    if (bc_pointer || this->_force_boundary_execution)
    {
      if (bc_pointer)
        bc_pointer->setupFaceData(&fi, this->_current_face_type);
      const auto rhs_contribution = this->computeBoundaryRHSContribution(*bc_pointer);

      const bool elem_side = this->_current_face_type == FaceInfo::VarFaceNeighbors::ELEM;
      const auto & info = elem_side ? *fi.elemInfo() : *fi.neighborInfo();
      const Real w = weight(elem_side ? fi.elemPtr() : fi.neighborPtr());
      const auto dof = info.dofIndices()[this->_sys_num][this->_var_num];
      (*this->_linear_system.rhs).add(dof, w * rhs_contribution);
    }
  }
}

using LinearFVTurbulentBlendedAdvection = LinearFVTurbulentBlended<LinearFVTurbulentAdvection>;
using LinearFVTurbulentBlendedDiffusion = LinearFVTurbulentBlended<LinearFVTurbulentDiffusion>;
using LinearFVTurbulentBlendedTimeDerivative =
    LinearFVTurbulentBlended<LinearFVTurbulentTimeDerivative>;
