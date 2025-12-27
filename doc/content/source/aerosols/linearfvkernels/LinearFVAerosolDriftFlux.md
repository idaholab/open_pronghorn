# LinearFVAerosolDriftFlux

## Overview

`LinearFVAerosolDriftFlux` adds the aerosol *drift-flux* contribution
$\rho \, \mathbf{U}_\text{drift} \cdot \mathbf{n} \, C_f$
to the linear finite-volume transport equation for a scalar quantity $C$
(e.g. aerosol mass or volume fraction).

The kernel is designed to be used with a `MooseLinearVariableFVReal` variable and
contributes both to the system matrix and, implicitly, to the right-hand side of
the linear system through face-based fluxes. It is typically combined with other
`LinearFVKernels` such as advection, diffusion, and source terms, and with
`LinearFV` aerosol boundary conditions (e.g. `LinearFVAerosolDepositionBC`).

The drift velocity $\mathbf{U}_\text{drift}$ is computed internally from:

- *Gravitational settling* (Stokes settling with Cunningham slip correction), and
- *Thermophoresis* (particle motion in the direction of $-\nabla T$).

All fluid properties are provided as *functors* (`rho`, `mu`, `T`), so this kernel
can be coupled to a wide variety of flow / heat-transfer models.

## Governing equations and physical model

For a transported aerosol scalar $C$ (for instance, a number, mass, or volume
fraction), the conservation equation in drift-flux form can be written as

\begin{equation}
\frac{\partial (\rho C)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} C)
+ \nabla \cdot \bigl(\rho \, \mathbf{U}_\text{drift} \, C\bigr)
= S,
\end{equation}

where

- $\rho$ is the carrier-gas density,
- $\mathbf{u}$ is the carrier-gas velocity (handled by other kernels), and
- $\mathbf{U}_\text{drift}$ is the relative drift velocity of particles with respect to the gas.

This kernel is responsible only for the drift-flux term

\begin{equation}
\nabla \cdot \bigl(\rho \, \mathbf{U}_\text{drift} \, C\bigr),
\end{equation}

while the advective term $\nabla \cdot (\rho \mathbf{u} C)$ and diffusive terms
are handled by separate `LinearFVKernels`.

### Gravitational settling

The gravitational drift component is computed assuming small particles in the
Stokes regime with a Cunningham slip correction. For spherical particles with

- diameter $d_p$ (`particle_diameter`),
- material density $\rho_p$ (`particle_density`),
- dynamic viscosity of the carrier gas $\mu$,
- gas mean free path $\lambda$ (`mean_free_path`), and
- gravitational acceleration vector $\mathbf{g}$ (`gravity`),

the Knudsen-number ratio and Cunningham slip factor are

\begin{align}
\text{Kn} &= \frac{\lambda}{d_p},\\
C_\mu &= 1 + \text{Kn} \left(2.34 + 1.05 \, e^{-0.39 d_p / \lambda}\right).
\end{align}

The particle relaxation time is

\begin{equation}
\tau_p = \frac{\rho_p d_p^2 C_\mu}{18 \, \mu},
\end{equation}

and the gravitational settling velocity is then

\begin{equation}
\mathbf{v}_g = \tau_p \, \mathbf{g}.
\end{equation}

In this kernel, $\mathbf{v}_g$ is evaluated at each face using the face-interpolated
viscosity $\mu_f$ and the user-specified particle properties and mean free path.

### Thermophoretic drift

The thermophoretic drift component moves particles down temperature gradients.
For a face temperature $T_f$, density $\rho_f$, viscosity $\mu_f$, and
temperature gradient $\nabla T_f$, the model implemented here is

\begin{equation}
\mathbf{v}_\text{th} = - k_\text{th} \, \frac{\mu_f}{\rho_f T_f} \, \nabla T_f,
\end{equation}

where $k_\text{th}$ is a thermophoretic coefficient.

The coefficient $k_\text{th}$ may be specified directly as a constant `k_th_const`, or,
if `use_constant_kth = false`, a simple diameter-based approximation is used
(currently a placeholder that can be refined as needed).

The temperature field and its gradient are obtained from the functor `T`.
The gradient is projected into the correct spatial dimension (1D, 2D, or 3D)
according to the problem mesh.

### Total drift velocity and face flux

The total drift velocity at a face is

\begin{equation}
\mathbf{U}_\text{drift} = \mathbf{v}_g + \mathbf{v}_\text{th}.
\end{equation}

The *drift mass flux* through a face with unit normal $\mathbf{n}$ is

\begin{equation}
m_{\text{drift}, f}
= \rho_f \, \mathbf{U}_\text{drift,f} \cdot \mathbf{n}_f,
\end{equation}

and the scalar flux is

\begin{equation}
F_f = m_{\text{drift}, f} \, C_f,
\end{equation}

where $C_f$ is the scalar value at the face, obtained using the selected
face interpolation scheme (see `advected_interp_method` below).

### Discretization and matrix contributions

`LinearFVAerosolDriftFlux` is a `LinearFVFluxKernel`, so it operates on mesh faces.
For each internal face, the face value $C_f$ is interpolated from the two
adjacent cell-centered values,

\begin{equation}
C_f = \alpha \, C_\text{elem} + (1 - \alpha) \, C_\text{neigh},
\end{equation}

where the interpolation weight $\alpha$ is chosen according to the
`advected_interp_method` (e.g. upwind, average, high-order schemes) and the sign
of the drift mass flux.

The resulting flux $F_f$ contributes linear terms to the residuals of the two
cells sharing the face. Because the flux is linear in $C$, this kernel only
assembles *matrix contributions*; the explicit right-hand-side contributions are
zero. This is reflected in the C++ implementation where

- `computeElemRightHandSideContribution()` and
- `computeNeighborRightHandSideContribution()`

both return zero.

### Boundary faces

By default, the parameter `force_boundary_execution = true` is used so that the
kernel executes on external boundaries even if no `LinearFVDirichletBC` is applied.
This allows the drift flux to contribute on open boundaries using the same
interpolation logic as for internal faces.

On boundaries where an explicit aerosol deposition boundary condition is applied
(e.g. `LinearFVAerosolDepositionBC`), the deposition BC provides an additional
normal flux, while `LinearFVAerosolDriftFlux` still supplies the bulk drift flux
on faces adjacent to the wall. Users should ensure that the combination of
kernels and boundary conditions reflects the intended physical model and that
flux terms are not double-counted.

!syntax parameters /LinearFVKernels/LinearFVAerosolDriftFlux

!syntax inputs /LinearFVKernels/LinearFVAerosolDriftFlux

!syntax children /LinearFVKernels/LinearFVAerosolDriftFlux
