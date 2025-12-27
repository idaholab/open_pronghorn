# LinearFVAerosolDepositionBC

## Overview

`LinearFVAerosolDepositionBC` is a finite-volume Robin boundary condition for aerosol
deposition at walls. It is intended for use with `MooseLinearVariableFVReal` variables
that satisfy an advection–diffusion(-drift) equation for an aerosol scalar $C$
(e.g. number, mass, or volume fraction).

This boundary condition:

- computes the *deposition velocity* $V_d$ using the Lai & Nazaroff (2000)
  semi-empirical model, based on friction velocity, Brownian diffusion, and
  gravitational settling, and
- enforces the Robin condition

\begin{equation}
\nabla C \cdot \mathbf{n} + \frac{V_d}{D} C = 0,
\end{equation}

so that the *diffusive flux into the wall* is

\begin{equation}
- D \, \nabla C \cdot \mathbf{n} = V_d \, C,
\end{equation}

where $D$ is the diffusion coefficient used in the bulk PDE.

Internally, `LinearFVAerosolDepositionBC` derives from
`LinearFVAdvectionDiffusionFunctorRobinBC` and overrides the Robin coefficients
$\alpha$, $\beta$, and $\gamma$ as

- $\alpha = 1$,
- $\beta = V_d / D$,
- $\gamma = 0$.

This BC is designed to be used together with drift-flux and diffusion kernels such as
`LinearFVAerosolDriftFlux` and `LinearFV{Anisotropic}Diffusion`, which model the
bulk aerosol transport.

## Robin form and implementation

### Robin form

The base class `LinearFVAdvectionDiffusionFunctorRobinBC` implements a canonical Robin
boundary condition for a scalar field $\phi$ of the form

\begin{equation}
\alpha \, \nabla \phi_b \cdot \hat{n}_b + \beta \, \phi_b = \gamma,
\end{equation}

where $\phi_b$, $\nabla \phi_b$, and $\hat{n}_b$ are the field value,
gradient, and outward normal at the boundary.

`LinearFVAerosolDepositionBC` specializes this for an aerosol scalar $C$ and the
Lai & Nazaroff deposition model by choosing

\begin{equation}
\alpha = 1, \qquad
\beta = \frac{V_d}{D}, \qquad
\gamma = 0.
\end{equation}

The resulting condition,

\begin{equation}
\nabla C \cdot \mathbf{n} + \frac{V_d}{D} C = 0,
\end{equation}

can be interpreted physically as a *sink* of aerosols at the wall, with a normal
diffusive flux proportional to the boundary concentration:

\begin{equation}
- D \, \nabla C \cdot \mathbf{n} = V_d \, C.
\end{equation}

### Coefficient evaluation

The Robin coefficients are provided through three virtual methods:

- `getAlpha(...)` $\rightarrow \alpha$,
- `getBeta(...)` $\rightarrow \beta$,
- `getGamma(...)` $\rightarrow \gamma$.

For `LinearFVAerosolDepositionBC`:

- `getAlpha` always returns `1.0`,
- `getGamma` always returns `0.0`,
- `getBeta` computes

  1. the *effective diffusion coefficient* $D$ at the boundary face using
     the functor `diffusion_coeff`, and

  2. the *deposition velocity* $V_d$ using
     `computeDepositionVelocity(...)`.

  If either $D \le 0$ or $V_d \le 0$, the method returns `0.0`
  (no deposition), which corresponds to a homogeneous Neumann condition for $C$.

The diffusion coefficient $D$ must be consistent with the coefficient used in the
bulk PDE, e.g. the sum of Brownian and turbulent diffusivities (`D_eff`) provided by
`AerosolDiffusivityFunctorMaterial` or a similar material.

## Deposition velocity model

The deposition velocity $V_d$ is computed following Lai & Nazaroff (2000) using
local wall quantities obtained from functors and the face normal.

### Input quantities at the wall

At each boundary face, the following quantities are evaluated from user-supplied
functors:

- $\rho$: carrier fluid density (`rho`),
- $\mu$: carrier fluid dynamic viscosity (`mu`),
- $T$: carrier fluid temperature (`T`),
- $u_*$: friction velocity at the wall (`u_star`).

In addition, the BC uses particle properties and gas mean free path from parameters:

- $d_p$: particle diameter (`particle_diameter`),
- $\rho_p$: particle density (`particle_density`),
- $\lambda$: mean free path (`mean_free_path`),
- $\mathbf{g}$: gravity vector (`gravity`).

If any of $\mu$, $\rho$, $T$, or $u_*$ are non-positive, the
computed deposition velocity is set to zero and the BC reduces to a zero-flux
condition for $C$.

### Gravitational settling

The gravitational settling velocity vector $\mathbf{v}_{gs}$ at the wall is
computed using the same Stokes–Cunningham formulation as in the drift-flux kernel.

The Cunningham slip correction factor is

\begin{equation}
C_c = 1 + \frac{\lambda}{d_p}
\left(2.34 + 1.05 \, e^{-0.39 d_p / \lambda}\right).
\end{equation}

The particle relaxation time is

\begin{equation}
\tau_p = \frac{\rho_p d_p^2 C_c}{18 \, \mu},
\end{equation}

and the gravitational settling velocity is

\begin{equation}
\mathbf{v}_{gs} = \tau_p \, \mathbf{g}.
\end{equation}

On the wall, only the *wall-normal component* of settling is relevant:

\begin{equation}
V_{gs} = \left|\mathbf{v}_{gs} \cdot \mathbf{n}\right|.
\end{equation}

The direction and magnitude of $\mathbf{g}$, combined with the face normal
$\mathbf{n}$, determine whether the face is treated as a floor, a ceiling,
or a vertical wall.

### Surface orientation: floor, ceiling, or vertical wall

Define

\begin{equation}
\text{alignment} = \frac{\mathbf{n} \cdot \mathbf{g}}{\lVert\mathbf{g}\rVert},
\end{equation}

and a cutoff $\text{horiz\_cutoff} \approx 0.5$ (corresponding to
about 60° from horizontal).

- If $|\text{alignment}| > \text{horiz\_cutoff}$, the surface is treated
  as *approximately horizontal*:
  - $\text{alignment} > 0$: wall normal roughly aligned with gravity
    $\Rightarrow$ *floor*,
  - $\text{alignment} < 0$: wall normal roughly opposite gravity
    $\Rightarrow$ *ceiling*.

- If $|\text{alignment}| \le \text{horiz\_cutoff}$, the surface is treated
  as *vertical*.

If the gravity magnitude is zero (or extremely small), the BC falls back to the
vertical-wall model for $V_d$.

### Horizontal surfaces: floor and ceiling

For horizontal surfaces, the Lai & Nazaroff model gives closed-form expressions
for the deposition velocity $V_d$ in terms of the wall-normal settling
velocity component $V_{gs}$ and the friction velocity $u_*$.

#### Floor (normal aligned with gravity)

For floors, with $\mathbf{n}$ aligned with $\mathbf{g}$, the model is

\begin{equation}
V_{d,f} = \frac{V_{gs}}{\exp(V_{gs}/u_*) - 1}.
\end{equation}

To avoid numerical issues:

- if $V_{gs}/u_*$ is very small, $V_{d,f}$ is limited to $u_*$, and
- if the denominator $\exp(V_{gs}/u_*) - 1$ is nearly zero, the BC also
  returns $V_{d,f} \approx u_*$.

#### Ceiling (normal opposite gravity)

For ceilings, with $\mathbf{n}$ approximately opposite $\mathbf{g}$,
the model is

\begin{equation}
V_{d,c} = \frac{V_{gs}}{1 - \exp(-V_{gs}/u_*)}.
\end{equation}

Analogous small-argument and near-singularity safeguards are applied, with
$V_{d,c}$ limited to $u_*$ as $V_{gs}/u_* \to 0$.

### Vertical walls

For vertical walls, Lai & Nazaroff propose an expression based on an integral
over the near-wall turbulent boundary layer:

\begin{equation}
V_{d,v} = \frac{u_*}{I},
\end{equation}

where the integral $I$ depends on particle size and flow properties through
the Schmidt number and a dimensionless wall coordinate $r^+$.

The implementation follows their formulation closely:

1. Compute the *Brownian diffusivity* (with slip correction):

   \begin{equation}
   D_B = \frac{k_B T C_c}{3 \pi \mu d_p},
   \end{equation}

   where $k_B$ is the Boltzmann constant and $C_c$ is the Cunningham
   slip factor defined above.

2. Compute the kinematic viscosity

   \begin{equation}
   \nu = \frac{\mu}{\rho},
   \end{equation}

   and the Schmidt number

   \begin{equation}
   Sc = \frac{\nu}{D_B}.
   \end{equation}

3. Compute the dimensionless radius

   \begin{equation}
   r^+ = \frac{d_p u_*}{2 \nu}.
   \end{equation}

4. Define $A_0 = 10.92\,Sc^{-1/3}$ and the functions

   - $a$: evaluated at $y^+ = 4.3$,
   - $b$: evaluated at $r^+$,

   with logarithmic and arctangent terms as given in Lai & Nazaroff. The
   implementation reproduces their table-based approximations.

5. Construct the integral

   \begin{equation}
   I = 3.64\,Sc^{2/3}(a - b) + 39.
   \end{equation}

If $I \le 0$ or not finite, $I$ is set to zero and the vertical-wall
deposition velocity $V_{d,v}$ is taken to be zero.

Finally,

\begin{equation}
V_{d,v} =
\begin{cases}
u_* / I, & I > 0, \\
0, & \text{otherwise.}
\end{cases}
\end{equation}

### Clipping and robustness

After computing $V_d$ as floor, ceiling, or vertical-wall velocity, the BC
enforces the following sanity checks:

- if $V_d$ is negative or not finite, it is set to zero;
- if diffusion coefficient $D \le 0$ at the boundary face, $\beta = 0$.

These checks ensure that the boundary condition is always well-posed and avoids
NaNs or non-physical negative deposition velocities.

!syntax parameters /LinearFVBCs/LinearFVAerosolDepositionBC

!syntax inputs /LinearFVBCs/LinearFVAerosolDepositionBC

!syntax children /LinearFVBCs/LinearFVAerosolDepositionBC
