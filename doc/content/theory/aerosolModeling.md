# Aerosol Modeling with Drift-Flux in OpenPronghorn

This page documents the aerosol-transport modeling capabilities implemented in
OpenPronghorn based on a *drift-flux* formulation. This aerosol transport formulation
is intended to be used with a baseline RANS models for determining the velocity field
and background turbulent diffusivity.

We focus on:

- the *drift-flux transport equation* for an aerosol scalar,
- the *gravitational settling* and *thermophoretic* drift velocities,
- the *Brownian* and *turbulent* diffusion models,
- the *wall deposition* model of Lai & Nazaroff,
- and how these are realized in OpenPronghorn via
  `LinearFVAerosolDriftFlux`, `AerosolDiffusivityFunctorMaterial`, and
  `LinearFVAerosolDepositionBC`.

For general background on aerosol dynamics, see standard references on aerosol science
(e.g., Hinds, *Aerosol Technology*) and indoor aerosol modeling
(e.g., [!cite](lai2000deposition), [!cite](nazaroff1993indoor). The drift-flux
implementation here follows [!cite](gray2025driftflux) for CFD modeling of
aerosol transport and deposition in ventilated domains.


## Introduction

OpenPronghorn applications frequently need to represent transport and deposition of
small, dilute particles (aerosols) in complex flows, for example:

- fission-product aerosols in reactor containment,
- smoke and soot transport,
- dust and particulate transport in ventilation systems,
- surrogate aerosols in experimental facilities.

Directly resolving the dispersed phase with Lagrangian particle tracking or two-fluid
Eulerian–Eulerian models can be expensive and intrusive when only *mean aerosol
concentration* and deposition patterns are needed. The *drift-flux* approach used
here is attractive because:

- it treats aerosols as a *passive scalar* (mass or volume fraction) with additional
  drift and diffusion terms;
- it leverages the existing *finite-volume Navier–Stokes* infrastructure and turbulence
  models;
- it is robust and comparatively inexpensive, making it suitable for large, complex
  reactor-relevant domains.

The key idea is to write a transport equation for a scalar aerosol quantity $C$ that
includes:

- advection by the carrier-gas velocity $\mathbf{u}$,
- *drift flux* due to gravitational settling and thermophoresis,
- *diffusion* due to Brownian motion and turbulent mixing,
- *wall deposition* via a physically based Robin boundary condition.


## Drift-flux formulation of aerosol transport

### Definition of the transported scalar

We introduce a scalar aerosol quantity $C$ defined at each point in the carrier gas.
Different choices are possible:

- *mass fraction* $Y_p = \rho_p^*/\rho$ of particulate material,
- *number density* $n$ of particles (per unit volume),
- *volume fraction* $\alpha_p$ of particles in the gas.

In OpenPronghorn, the aerosol modeling infrastructure is agnostic to which specific
scalar is used, provided the source/sink terms and post-processing are consistent.
For the purposes of this document, we refer to $C$ generically as "aerosol
concentration".

We assume:

- the aerosol phase is *dilute*, so its effect on the carrier-gas density and
  viscosity is negligible;
- all particles are *monodisperse*, with a single diameter $d_p$ and density
  $\rho_p$;
- particles are *small and heavy* enough that their dynamics can be approximated by
  a relaxation-time model with *Stokes drag* and a *Cunningham* slip correction;
- inter-particle interactions (coagulation, breakup) and thermodynamic phase change
  are neglected in the base model.

### Drift-flux transport equation

Let $\rho$ denote the carrier-gas density and $\mathbf{u}$ its mean velocity,
as determined by the Navier–Stokes solution (possibly RANS). The drift-flux transport
equation for $C$ is written as

\begin{equation}
\frac{\partial (\rho C)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} C)
+ \nabla \cdot \bigl(\rho \mathbf{U}_\text{drift} C\bigr)
- \nabla \cdot \bigl(\rho D \nabla C\bigr)
= S_C,
\label{eq:aerosol-transport}
\end{equation}

where

- $\rho \mathbf{u} C$ is *bulk advection* by the gas velocity,
- $\rho \mathbf{U}_\text{drift} C$ is the *drift flux* relative to the gas,
- $\rho D \nabla C$ is the diffusive flux,
- $S_C$ represents sources and sinks (e.g., injection, bulk removal, reactions).

The *drift velocity* $\mathbf{U}_\text{drift}$ is modeled as the sum of

\begin{equation}
\mathbf{U}_\text{drift} = \mathbf{v}_g + \mathbf{v}_\text{th},
\end{equation}

where

- $\mathbf{v}_g$ is the *gravitational settling velocity*,
- $\mathbf{v}_\text{th}$ is the *thermophoretic drift velocity*.

The *effective diffusivity* $D$ is modeled as

\begin{equation}
D = D_B + D_t,
\end{equation}

the sum of Brownian diffusivity $D_B$ and turbulent diffusivity $D_t$.

In OpenPronghorn, Eq. \eqref{eq:aerosol-transport} is discretized with finite volumes
using:

- `LinearFVAdvection` for the $\rho \mathbf{u} C$ term,
- `LinearFVAerosolDriftFlux` for the $\rho \mathbf{U}_\text{drift} C$ term,
- `LinearFVDiffusion` (or `LinearFVAnisotropicDiffusion`) with diffusion coefficient
  from `AerosolDiffusivityFunctorMaterial` for $\rho D \nabla C$,
- optional source/sink kernels for $S_C$.


## Gravitational settling model

### Stokes drag and relaxation time

We assume the particles are small enough that their motion relative to the gas can be
described by Stokes drag, corrected for slip effects. The *particle relaxation time*
$\tau_p$ is

\begin{equation}
\tau_p = \frac{\rho_p d_p^2 C_c}{18 \, \mu},
\label{eq:taup}
\end{equation}

where

- $d_p$ is the particle diameter,
- $\rho_p$ is the particle density,
- $\mu$ is the gas dynamic viscosity,
- $C_c$ is the *Cunningham slip correction*.

The Cunningham correction accounts for the fact that as $d_p$ becomes comparable to
the gas mean free path $\lambda$, the no-slip boundary condition at the particle
surface becomes less accurate. In the implementation we use

\begin{equation}
C_c = 1 + \frac{\lambda}{d_p}
\left( 2.34 + 1.05 \, e^{-0.39 d_p / \lambda} \right),
\label{eq:cunningham}
\end{equation}

where $\lambda$ is a user-specified *mean free path* (`mean_free_path` parameter).
For typical air at standard conditions, $\lambda \approx 0.066 \,\mu\text{m}$.

### Settling velocity

Given $\tau_p$, the *gravitational settling velocity* is

\begin{equation}
\mathbf{v}_g = \tau_p \, \mathbf{g},
\label{eq:vg}
\end{equation}

where $\mathbf{g}$ is the gravitational acceleration vector. The sign and magnitude
of $\mathbf{g}$ determine the direction of settling; in most applications
$\mathbf{g} = (0,0,-9.81)\,\text{m/s}^2$.

Equation \eqref{eq:vg} is evaluated locally at each cell face using interpolated values
of $\mu$ and the user-specified particle properties. In OpenPronghorn this is
implemented in `LinearFVAerosolDriftFlux` using functors `mu`, `rho`, and the parameters
`particle_diameter`, `particle_density`, and `mean_free_path`.


## Thermophoretic drift model

### Physical mechanism

Thermophoresis is the migration of particles in the presence of a temperature gradient.
Gas molecules from the hotter side have higher thermal velocities and impart larger
momentum to the particle surface than those from the colder side, resulting in a net
force pushing particles from hot to cold regions.

Thermophoresis is particularly important near *strongly heated or cooled walls*,
where $\nabla T$ can be large and cause significant drift toward or away from the
surface, depending on the sign convention.

### Drift velocity expression

In the implementation, the thermophoretic drift velocity is modeled as

\begin{equation}
\mathbf{v}_\text{th} = - k_\text{th}
\, \frac{\mu}{\rho T} \, \nabla T,
\label{eq:vth}
\end{equation}

where

- $k_\text{th}$ is a *thermophoretic coefficient*,
- $T$ is the local gas temperature,
- $\mu$ and $\rho$ are the gas viscosity and density,
- $\nabla T$ is the local temperature gradient.

The form in Eq. \eqref{eq:vth} is consistent with the drift-flux model in
Gray et al. [!cite](gray2025driftflux), which writes thermophoretic drift in terms of
$\nu / (\rho T)$ and $\nabla T$, where $\nu = \mu / \rho$ is the kinematic
viscosity.

### Thermophoretic coefficient

The coefficient $k_\text{th}$ depends on gas properties, particle properties, and
Knudsen number. In many engineering applications it is approximated as a *constant*
for a given particle size and carrier gas. OpenPronghorn therefore supports:

- a *constant* value `k_th_const` (default $k_\text{th} = 0.5$),
- a placeholder for future diameter- or Knudsen-dependent formulas.

When `use_constant_kth = true`, Eq. \eqref{eq:vth} is evaluated with
$k_\text{th} = \texttt{k\_th\_const}$. The scalar coefficient and the gradient
$\nabla T$ are evaluated at each face using the `T` functor.

Thermophoresis and gravitational settling are simply added:

\begin{equation}
\mathbf{U}_\text{drift} = \mathbf{v}_g + \mathbf{v}_\text{th}.
\end{equation}


## Diffusion: Brownian and turbulent

### Brownian diffusion $D_B$

Small aerosol particles undergo Brownian motion due to collisions with gas molecules.
Assuming Stokes drag with slip correction, a Stokes–Einstein-type expression for the
Brownian diffusivity is

\begin{equation}
D_B = \frac{k_B T C_c}{3 \pi \mu d_p},
\label{eq:DB}
\end{equation}

where

- $k_B$ is the Boltzmann constant,
- $T$, $\mu$ and $d_p$ are as defined above,
- $C_c$ is the Cunningham correction given by Eq. \eqref{eq:cunningham}.

This expression increases with temperature and decreases with particle size and
viscosity. For particles larger than a few microns, Brownian diffusion is typically
small compared to turbulent diffusion, but it remains important for submicron aerosols.

In `AerosolDiffusivityFunctorMaterial`, Eq. \eqref{eq:DB} is evaluated using the
functors `T` and `mu`, the material parameter `particle_diameter`, and the mean free
path `mean_free_path`. A small positive lower bound $D_{B,\min}$ is applied to avoid
exactly zero diffusivity when inputs become non-physical (e.g., due to numerical noise).

### Turbulent diffusion $D_t$

Turbulence enhances mixing at scales larger than the molecular mean free path.
We model turbulent diffusion using an *eddy-diffusivity* approximation based on the
turbulent viscosity $\mu_t$:

\begin{equation}
D_t = \frac{\mu_t}{\rho \, Sc_t},
\label{eq:Dt}
\end{equation}

where

- $\mu_t$ is the turbulent viscosity (e.g., from a $k$–$\epsilon$ model),
- $\rho$ is the gas density,
- $Sc_t$ is the *turbulent Schmidt number*, typically of order unity.

In OpenPronghorn, $\mu_t$ is provided by the turbulence model via a functor
(`NS::mu_t`), and `AerosolDiffusivityFunctorMaterial` computes $D_t$ using Eq.
\eqref{eq:Dt}. A small lower bound $D_{t,\min}$ is again used to maintain numerical
stability when $\mu_t$ or $\rho$ approach zero.

The parameter `Sc_t` is user-configurable. The default value (1.0) is commonly used for
neutral, shear-dominated flows. Lower values (e.g., 0.7) increase turbulent diffusion, and
higher values decrease it.

### Effective diffusivity $D_\text{eff}$

The *effective* aerosol diffusivity is defined as the sum

\begin{equation}
D_\text{eff} = D_B + D_t.
\label{eq:Deff}
\end{equation}

In the transport equation \eqref{eq:aerosol-transport}, it is this effective diffusivity
that appears in the diffusive flux term. `AerosolDiffusivityFunctorMaterial` exposes
three functor properties:

- `D_B` for Brownian diffusion,
- `D_t` for turbulent diffusion,
- `D_eff` for the total effective diffusion, Eq. \eqref{eq:Deff}.

Kernels and boundary conditions can then select the appropriate functor for their
coefficients. For example, `LinearFVDiffusion` and `LinearFVAerosolDepositionBC`
typically both use `D_eff` to ensure consistency between bulk diffusion and wall
deposition.


## Wall deposition modeling

### Physical picture

Aerosol deposition at solid surfaces results from several mechanisms:

- *Brownian diffusion* bringing particles into contact with the wall,
- *turbulent diffusion* transporting particles across the boundary layer,
- *gravitational settling* toward horizontal surfaces (floors, ceilings),
- *thermophoresis* driving particles toward hot or cold walls,
- possibly other mechanisms (electrophoresis, interception, impaction) not included
  in the base model.

Instead of explicitly resolving all near-wall processes, we use a *deposition
velocity* model: the net normal particle flux at the wall is proportional to the
local concentration,

\begin{equation}
J_\text{dep} = V_d \, C_w,
\end{equation}

where $C_w$ is the aerosol concentration at the wall and $V_d$ is the *deposition
velocity* (units of m/s). This is implemented as a Robin boundary condition on $C$.

### Lai & Nazaroff model

OpenPronghorn adopts the semi-empirical model of Lai & Nazaroff [!cite](lai2000deposition)
for deposition in turbulent boundary layers. The model expresses $V_d$ in terms of

- friction velocity $u_*$,
- particle size $d_p$,
- Brownian diffusivity $D_B$,
- kinematic viscosity $\nu = \mu/\rho$,
- gravitational settling velocity $\mathbf{v}_{gs}$,
- and surface orientation relative to gravity.

The normal aerosol flux at the wall is modeled as

\begin{equation}
- D \, \nabla C \cdot \mathbf{n} = V_d \, C,
\end{equation}

with $D$ taken as the same effective diffusivity $D_\text{eff}$ used in the bulk
transport equation. In Robin form this becomes

\begin{equation}
\nabla C \cdot \mathbf{n} + \frac{V_d}{D} C = 0.
\end{equation}

`LinearFVAerosolDepositionBC` enforces this condition by setting the Robin coefficients
$\alpha = 1$, $\beta = V_d/D$, and $\gamma = 0$.

### Surface orientation and gravitational settling

The model distinguishes *floor*, *ceiling*, and *vertical* surfaces based on the
alignment between the wall-normal vector $\mathbf{n}$ and gravity $\mathbf{g}$. Define

\begin{equation}
\text{alignment} = \frac{\mathbf{n} \cdot \mathbf{g}}{\|\mathbf{g}\|}.
\end{equation}

Using a cutoff $|\text{alignment}| > \text{horiz\_cutoff} \approx 0.5$, we classify:

- *Horizontal surfaces*:
  - $\text{alignment} > 0$: $\mathbf{n}$ aligned with $\mathbf{g}$ – *floor*,
  - $\text{alignment} < 0$: $\mathbf{n}$ opposite $\mathbf{g}$ – *ceiling*.
- *Vertical surfaces*:
  - $|\text{alignment}| \le \text{horiz\_cutoff}$.

The *wall-normal component* of gravitational settling is

\begin{equation}
V_{gs} = |\mathbf{v}_{gs} \cdot \mathbf{n}|,
\end{equation}

with $\mathbf{v}_{gs}$ given by Eqs. \eqref{eq:taup}–\eqref{eq:vg}.

When $\|\mathbf{g}\|$ is effectively zero, gravitational settling is neglected and
all surfaces are treated as "vertical" in the sense of the model.

### Floors and ceilings

For essentially horizontal surfaces, Lai & Nazaroff provide closed-form expressions for
$V_d$ in terms of $V_{gs}$ and friction velocity $u_*$.

- *Floor* (normal aligned with gravity):

  \begin{equation}
  V_{d,f} = \frac{V_{gs}}{\exp(V_{gs}/u_*) - 1}.
  \end{equation}

- *Ceiling* (normal opposite gravity):

  \begin{equation}
  V_{d,c} = \frac{V_{gs}}{1 - \exp(-V_{gs}/u_*)}.
  \end{equation}

In the implementation, small-argument limits are handled carefully to avoid division by
nearly zero denominators; when $V_{gs}/u_*$ is very small, both expressions tend
smoothly to $V_d \approx u_*$.

### Vertical walls

For *vertical* surfaces, Lai & Nazaroff derive an expression based on integrating the
balance of diffusion and turbulent transport across the near-wall region. The deposition
velocity is written as

\begin{equation}
V_{d,v} = \frac{u_*}{I},
\label{eq:Vd-vertical}
\end{equation}

where $I$ is an integral that depends on:

- the Schmidt number $Sc = \nu/D_B$,
- the dimensionless particle radius $r^+ = d_p u_* / (2 \nu)$,
- empirical functions $a$ and $b$ evaluated at specific wall coordinates.

The implementation follows Lai & Nazaroff's expressions:

1. Compute Brownian diffusivity $D_B$ and $\nu = \mu/\rho$.
2. Form the Schmidt number $Sc = \nu/D_B$.
3. Compute $r^+ = d_p u_* / (2 \nu)$.
4. Evaluate auxiliary quantities (e.g., $A_0 = 10.92\,Sc^{-1/3}$, functions $a$
   and $b$ at $y^+ = 4.3$ and $r^+$).
5. Assemble the integral
   \begin{equation}
   I = 3.64\,Sc^{2/3}(a - b) + 39.
   \end{equation}

If $I$ is positive and finite, Eq. \eqref{eq:Vd-vertical} is used; otherwise, the
deposition velocity is set to zero. This provides a robust expression for deposition on
vertical walls that smoothly connects regimes dominated by Brownian diffusion, turbulent
diffusion, and gravitational settling along the wall.

### Numerical robustness and clipping

After computing $V_d$ for floors, ceilings, or vertical walls, the implementation
enforces:

- $V_d \ge 0$,
- $V_d = 0$ if any required inputs ($u_*$, $\rho$, $\mu$, $T$) are non-positive
  or unphysical,
- $\beta = V_d/D = 0$ if $D \le 0$.

These checks guarantee a well-posed Robin boundary condition and avoid numerical
issues from negative or NaN coefficients.


## Numerical implementation in OpenPronghorn

### Building blocks and coupling

The aerosol drift-flux model is implemented as a set of modular MOOSE objects:

- `AerosolDiffusivityFunctorMaterial`
  computes and exposes `D_B`, `D_t`, and `D_eff` as functor properties, based on
  `rho`, `mu`, `T`, `mu_t`, `particle_diameter`, `mean_free_path`, and `Sc_t`.

- `LinearFVAerosolDriftFlux`
  adds the drift-flux term $\nabla \cdot (\rho \mathbf{U}_\text{drift} C)$ to the
  aerosol transport equation, with
  $\mathbf{U}_\text{drift} = \mathbf{v}_g + \mathbf{v}_\text{th}$ modeled as described
  above.

- `LinearFVAerosolDepositionBC`
  implements the Lai & Nazaroff deposition boundary condition
  $\nabla C \cdot \mathbf{n} + (V_d/D) C = 0$ using `diffusion_coeff` (usually `D_eff`)
  and functors for `rho`, `mu`, `T`, `u_star`.

These are combined with existing Navier–Stokes and turbulence objects:

- mean velocity $\mathbf{u}$ and pressure $p$ from the NS system,
- density, viscosity, and temperature functors (`NS::density`, `NS::mu`, `NS::T_fluid`),
- turbulent viscosity `NS::mu_t` from a RANS model (e.g. `kEpsilonViscosity`).

The aerosol concentration variable itself is typically declared as

```text
[Variables]
  [C_a]
    type = MooseLinearVariableFVReal
    solver_sys = 'aerosol_sys'
    initial_condition = 0.0
  []
[]
```

and all aerosol kernels and BCs refer to this variable.


### Discretization and flux structure

The finite-volume discretization in OpenPronghorn is face-based. For each internal
face between two cells, the total flux of $C$ can be written schematically as

\begin{equation}
F_f = F_f^\text{adv} + F_f^\text{drift} + F_f^\text{diff} + F_f^\text{src},
\end{equation}

where

- $F_f^\text{adv}$ is advection by $\mathbf{u}$ (`LinearFVAdvection`),
- $F_f^\text{drift}$ is drift flux by $\mathbf{U}_\text{drift}$
  (`LinearFVAerosolDriftFlux`),
- $F_f^\text{diff}$ is diffusive flux based on $D_\text{eff}$
  (`LinearFVDiffusion`),
- $F_f^\text{src}$ collects any additional user-defined contributions.

Each flux is constructed in a conservative fashion using interpolation of face values
and normals, and all contributions assemble into the same linear system for `C_a`.

On boundaries, `LinearFVAerosolDepositionBC` provides a flux corresponding to
$V_d C_w$. Other types of boundaries (inlets, outlets) can be specified using the
standard finite-volume BCs (`LinearFVDirichletBC`, `LinearFVNeumannBC`, etc.) and
additional custom Robin conditions if needed.


## Practical guidance for aerosol modeling

### When to use the drift-flux model

The drift-flux approach described here is appropriate when:

- aerosol volume fraction is small (dilute limit),
- particle–particle interactions and inertia are moderate,
- interest is primarily in *mean concentration fields* and *deposition patterns*,
- the gas-phase flow is resolved with RANS or laminar NS.

It may be less appropriate when:

- the particle size distribution is broad and strongly size-dependent,
- ballistic or inertial impaction dominates (very large Stokes number),
- coagulation, condensation, or thermal decomposition are important,
- strongly transient puff-like release behavior must be captured in detail.

In such cases, more advanced models (polydisperse multi-group drift-flux, Lagrangian
tracking, or fully coupled two-fluid models) may be required.


### Coupling with turbulence models

Because turbulent diffusion is proportional to $\mu_t$ via Eq. \eqref{eq:Dt},
the choice of turbulence model directly affects aerosol dispersion. Practical advice:

- For *simple internal flows* with moderate Reynolds numbers and weak separation,
  high-Re $k$–$\epsilon$ variants with wall functions are often sufficient.
- For *thermal–hydraulic applications with strong near-wall heat transfer*, low-Re
  or two-layer $k$–$\epsilon$ models can improve predictions of near-wall
  temperature gradients and velocity profiles, which in turn influence Brownian and
  thermophoretic transport.
- For *swirling or rotating flows*, realizable or curvature-corrected turbulence
  models help obtain better estimates of $\mu_t$ and mixing.

Always ensure that:

- the turbulence model and aerosol model use *consistent wall-distance fields* and
  *boundary lists* (e.g., for `WallDistanceAux` and wall BCs),
- the *friction velocity* $u_*$ used in `LinearFVAerosolDepositionBC` is consistent
  with the velocity and turbulence solution (e.g., obtained from wall shear stress or
  wall-function relations).


### Typical parameter ranges

Some practical ranges and considerations:

- *Particle diameter $d_p$*:
  The model is most appropriate for diameters from $\mathcal{O}(0.01\,\mu\text{m})$
  to a few $\mu\text{m}$. At very large sizes (tens of microns), gravitational
  settling becomes highly dominant and inertial effects may require more detailed
  modeling; at very small sizes (nanoparticles), Brownian diffusion dominates and
  accurate gas-property modeling becomes important.

- *Mean free path $\lambda$*:
  For air-like conditions, $\lambda \sim 0.05$–$0.1\,\mu\text{m}$. Sensitivity is
  strongest for $d_p \lesssim \mathcal{O}(\lambda)$; for larger particles the
  Cunningham factor approaches 1.

- *Turbulent Schmidt number $Sc_t$*:
  Values around 0.7–1.0 are commonly used for passive scalars in shear flows. Smaller
  values increase aerosol diffusion, larger values reduce it. If in doubt, start with
  `Sc_t = 1.0` and calibrate as needed against data.

- *Thermophoretic coefficient $k_\text{th}$*:
  Typical constant values range from 0.3–0.7 depending on regime. The default
  $k_\text{th} = 0.5$ is a reasonable choice for many aerosol–air systems.


### Mesh and time-step considerations

Accurate aerosol modeling requires sufficient resolution to capture:

- *major flow features* (jets, recirculation zones, shear layers),
- *temperature gradients*, especially near heated or cooled walls,
- near-wall regions where deposition occurs.

Specific guidance:

- Use a mesh that adequately resolves the *boundary layers* in which
  $u_*$ and near-wall gradients are determined; this typically follows the same
  criteria as for the chosen turbulence model (e.g., $y^+$ ranges).
- Ensure that the mesh is fine enough to resolve *thermal gradients* where
  thermophoresis is significant (e.g., near hot walls or cooled surfaces).
- Time-step constraints are usually driven by the *carrier flow*. The aerosol scalar
  tends to be less restrictive, but extreme drift velocities (e.g., large $d_p$ and
  strong gravity) may require smaller time steps for stability and accuracy in
  transient simulations.


### Validation and calibration

As with turbulence modeling, it is important to validate aerosol simulations against
experimental data or higher-fidelity simulations:

- For deposition, compare predicted *deposition velocities* or *wall loading
  profiles* against controlled experiments (e.g., channel or chamber tests).
- For bulk transport, compare predicted *concentration fields* or *breakthrough
  curves* against tracer measurements.
- Use these comparisons to calibrate parameters such as `Sc_t`, `k_th_const`, and
  any additional empirical coefficients introduced in extended models.

Document any changes to default parameters and the validation data used to justify them,
especially for safety-related applications.

## Example applications and typical usage

The table below summarizes several common application types, how to configure the
drift-flux aerosol model in OpenPronghorn, and typical outcomes of interest.

| Example application | How to set up the model | Outcomes of common interest |
| --- | --- | --- |
| Fission-product aerosols in reactor containment | Solve RANS (e.g. $k$-$\epsilon$) for gas flow and temperature; add a MooseLinearVariableFVReal aerosol scalar (for example C_a); use AerosolDiffusivityFunctorMaterial with air-like properties, appropriate particle_diameter and Sc_t; add LinearFVAdvection, LinearFVDiffusion with D_eff, and LinearFVAerosolDriftFlux; apply LinearFVAerosolDepositionBC on containment walls using u_star from the turbulence model. | Spatial and temporal evolution of aerosol concentration in the containment; wall deposition rates and loads on key structures; airborne versus deposited inventory for source-term and risk assessments. |
| Ventilation duct / HVAC aerosol transport | Solve steady-state RANS in ductwork (including bends, junctions, filters); define an aerosol scalar with inlet concentration or mass-flux boundary condition; use AerosolDiffusivityFunctorMaterial for Brownian and turbulent diffusion; add LinearFVAerosolDriftFlux (gravity may be important in horizontal runs); apply LinearFVAerosolDepositionBC on duct walls and special boundary conditions on filters as needed. | Outlet aerosol concentration and removal efficiency; deposition patterns along duct walls (hot spots, bends, vertical versus horizontal legs); time-integrated mass loading on filters and downstream components. |
| Smoke or soot dispersion in a ventilated room | Solve RANS with buoyancy and heat transfer for a room or multi-room domain; add an aerosol scalar with volumetric or inlet sources to represent smoke or soot; use AerosolDiffusivityFunctorMaterial with soot-relevant particle_diameter and air properties; include LinearFVAerosolDriftFlux (gravity drift and thermophoresis near hot surfaces); apply LinearFVAerosolDepositionBC on walls, floor, and ceiling using u_star from wall models. | Time evolution of smoke concentration in occupied regions; visibility-related metrics derived from soot concentration; deposition patterns on walls and ceilings for cleaning and inspection planning. |
| Bench-scale validation in a canonical geometry (channel or cavity) | Set up a simple 2D or 3D channel or cavity matching experimental conditions; run laminar or RANS Navier–Stokes with prescribed inlet profiles and thermal boundary conditions; configure AerosolDiffusivityFunctorMaterial and LinearFVAerosolDriftFlux with experimental particle_diameter, particle_density, and mean_free_path; apply LinearFVAerosolDepositionBC on selected walls with u_star consistent with measurements; prescribe inlet aerosol concentration as in the experiment. | Comparison of predicted versus measured concentration profiles; comparison of effective deposition velocities on walls; sensitivity of results to Sc_t, k_th_const, and particle size. |
| Dust and particulate transport in reactor cavities or tunnels | Solve RANS in a multi-region domain with recirculating flow and internal structures; add aerosol sources at leak or resuspension locations; use AerosolDiffusivityFunctorMaterial with dust-relevant properties and calibrated Sc_t; add LinearFVAerosolDriftFlux to capture settling in low-velocity regions and optionally thermophoresis with strong thermal gradients; apply LinearFVAerosolDepositionBC on concrete and steel surfaces. | Identification of stagnation and accumulation regions (dust hot spots); wall loading and deposition maps for maintenance and radiological protection; time scales for transport from leak locations to monitored or filtered regions. |


## Summary

The drift-flux aerosol model in OpenPronghorn provides a flexible, modular framework for
predicting aerosol transport and deposition in complex thermal–hydraulic systems. It
combines:

- a conservative finite-volume transport equation for an aerosol scalar,
- gravitational and thermophoretic drift velocities based on particle and gas properties,
- Brownian and turbulent diffusion via `AerosolDiffusivityFunctorMaterial`,
- wall deposition via the Lai & Nazaroff model in `LinearFVAerosolDepositionBC`,
- seamless coupling with existing Navier–Stokes and turbulence models.

By following the guidelines in this document—choosing appropriate particle properties,
turbulence models, boundary conditions, and numerical settings—users can simulate a wide
range of aerosol problems in reactor and thermal–hydraulic applications within the
OpenPronghorn framework.
