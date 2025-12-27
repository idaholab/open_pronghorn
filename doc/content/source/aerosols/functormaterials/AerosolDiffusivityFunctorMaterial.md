# AerosolDiffusivityFunctorMaterial

## Overview

`AerosolDiffusivityFunctorMaterial` is a `FunctorMaterial` that provides
aerosol diffusion coefficients as *functor properties* for use in finite-volume
and finite-element kernels.

It computes and exposes three scalar functors:

- `D_B` — Brownian (molecular) aerosol diffusivity,
- `D_t` — turbulent aerosol diffusivity,
- `D_eff` — effective diffusivity, defined as
  $D_\text{eff} = D_B + D_t$.

These functors can be used, for example, by

- `LinearFVDiffusion` / `LinearFVAnisotropicDiffusion` kernels to model
  aerosol diffusion in the bulk, and
- `LinearFVAerosolDepositionBC` (via its `diffusion_coeff` parameter) to
  supply a physically consistent diffusion coefficient in the wall deposition
  boundary condition.

All quantities are computed using local thermophysical properties obtained
from user-specified functors (`rho`, `mu`, `T`, `mu_t`) together with particle
size and gas mean free path.

## Physics and constitutive model

### Brownian diffusion `D_B`

Brownian diffusion represents the random motion of small aerosol particles
due to collisions with gas molecules. The material computes `D_B` using a
Stokes–Einstein relation corrected by the Cunningham slip factor,

\begin{equation}
D_B \propto \frac{k_B T \, C_c}{\mu \, d_p},
\end{equation}

where

- $k_B$ is the Boltzmann constant,
- $T$ is the local gas temperature,
- $\mu$ is the local dynamic viscosity,
- $d_p$ is the particle diameter (`particle_diameter`),
- $C_c$ is the Cunningham slip correction, a function of $d_p$ and
  the gas mean free path $\lambda$ (`mean_free_path`).

The Cunningham correction factor takes the form

\begin{equation}
C_c = 1 + \frac{\lambda}{d_p}
\left(2.34 + 1.05 \, e^{-0.39 d_p / \lambda}\right),
\end{equation}

which increases the diffusivity for small particles as their size approaches
the gas mean free path.

To ensure numerical robustness, a small positive lower bound
$D_{B,\min} > 0$ is applied so that `D_B` never becomes exactly zero.
If any of the required inputs (e.g. $T$, $\mu$, or $d_p$) are
non-physical (non-positive), the material falls back to this minimum value.

### Turbulent diffusion `D_t`

Turbulent diffusion represents enhanced mixing caused by turbulent eddies.
The material models `D_t` using an eddy-diffusivity approximation,

\begin{equation}
D_t = \frac{\mu_t}{\rho \; Sc_t},
\end{equation}

where

- $\mu_t$ is the turbulent (eddy) viscosity, provided via the functor
  `mu_t`,
- $\rho$ is the local fluid density, provided via the functor `rho`,
- $Sc_t$ is the turbulent Schmidt number (`Sc_t`).

If $\rho \le 0$ or $Sc_t \le 0$, the model deactivates turbulent
diffusion and returns a small positive lower bound $D_{t,\min}$ to keep
the coefficient non-zero and numerically stable.

Choosing $Sc_t \approx 1$ is common in near-neutral flows, but other
values can be set depending on the turbulence model and calibration data.

### Effective diffusivity `D_eff`

The effective aerosol diffusivity is defined simply as

\begin{equation}
D_\text{eff} = D_B + D_t.
\end{equation}

This functor property can be used directly in bulk diffusion kernels and
boundary conditions, ensuring that molecular and turbulent contributions are
treated consistently across the domain.

Internally, `D_eff` is computed by summing the `D_B` and `D_t` functor
properties provided by the same material, so no additional fluid or particle
parameters are required.

## Input parameters

`AerosolDiffusivityFunctorMaterial` extends the standard `FunctorMaterial`
interface. The most relevant parameters are summarized below.

### Required parameters

- `rho` (`MooseFunctorName`)
  Functor providing the carrier-gas density $\rho$. This is
  used in the turbulent diffusion model.

- `mu` (`MooseFunctorName`)
  Functor providing the carrier-gas dynamic viscosity $\mu$. This is
  used in the Brownian diffusion model and in forming the kinematic viscosity.

- `T` (`MooseFunctorName`)
  Functor providing the carrier-gas temperature $T$. This is used in the
  Brownian diffusion model.

- `mu_t` (`MooseFunctorName`)
  Functor providing the turbulent (eddy) viscosity $\mu_t$. This is
  used in the turbulent diffusion model. If you do not wish to include turbulent
  diffusion, you can supply a functor that returns zero or choose a very large
  `Sc_t` to effectively suppress $D_t$.

- `particle_diameter` (`Real`)
  Particle diameter $d_p$ used in both the Brownian model and the
  Cunningham slip correction. The material assumes monodisperse particles of
  this size.

### Optional parameters

- `mean_free_path` (`Real`)
  Gas mean free path $\lambda$ used in the Cunningham correction for
  Brownian diffusion. For air at standard conditions, $\lambda$ is typically
  $\mathcal{O}(10^{-7})$ m (tens of nanometers). The default value in the
  input syntax table reflects a typical choice for air.

- `Sc_t` (`Real`)
  Turbulent Schmidt number $Sc_t$ used in the turbulent diffusion model
  $D_t = \mu_t / (\rho Sc_t)$. The default is typically of order unity;
  see the syntax table below for the exact default value.

All other standard `FunctorMaterial` parameters (e.g. `block`, `execute_on`)
are also available and control where and when the material is evaluated.

## Provided functor properties

This material defines the following functor properties, which can be accessed by
other objects (kernels, BCs, materials) via their `MooseFunctorName` parameters:

- `D_B` — Brownian diffusivity of the aerosol,
- `D_t` — turbulent diffusivity of the aerosol,
- `D_eff` — effective diffusivity (Brownian + turbulent).

For example, in another object you can refer to the effective diffusivity using

```text
diffusion_coeff = D_eff
```

provided that `AerosolDiffusivityFunctorMaterial` is active on the same blocks.

!syntax parameters /FunctorMaterials/AerosolDiffusivityFunctorMaterial

!syntax inputs /FunctorMaterials/AerosolDiffusivityFunctorMaterial

!syntax children /FunctorMaterials/AerosolDiffusivityFunctorMaterial
