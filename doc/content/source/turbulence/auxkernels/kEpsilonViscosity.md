# kEpsilonViscosity

`kEpsilonViscosity` is an auxiliary kernel that computes the *turbulent dynamic viscosity*
$\mu_t$ used in the k–$\epsilon$ family of turbulence models.
It forms the closure term for the Reynolds stresses in the momentum equations:

\begin{equation}
\mu_t = \rho\, C_\mu\, k\, T,
\end{equation}

where the *turbulent time scale* $T$ and the *effective coefficient* $C_\mu$ depend on the
selected k–$\epsilon$ variant and the local flow state.

This object implements:

- Standard, Low-Re, Two-Layer, and Realizable viscosity formulations,
- Several wall treatments (Newton, incremental, linearized, non-equilibrium),
- Two-layer near-wall blending methods ([!cite](wolfstein1969velocity), [!cite](norris1975oneequation)),
- Low-Re damping functions for the StandardLowRe model,
- Realizable variable $C_\mu$,
- Optional scale limiting for the time scale,
- Bulk and near-wall formulations that match Menter-type corrections [!cite](menter1994two).


## Turbulent viscosity in the bulk region

Away from walls (or when `[bulk_wall_treatment=false]`), the turbulent viscosity is computed as:

\begin{equation}
\mu_t = \rho\, C_\mu^{\mathrm{eff}}\, k\, T,
\end{equation}

where:

- $k$ is the turbulent kinetic energy,
- $\rho$ is the fluid density,
- $C_\mu^{\mathrm{eff}}$ depends on the k–$\epsilon$ variant,
- $T$ is the turbulent time scale (possibly limited).

### Turbulent time scale

The base time scale is:

\begin{equation}
T_e = \frac{k}{\epsilon}.
\end{equation}

If `[scale_limiter=standard]`, the time scale is limited using:

\begin{equation}
T = \max\left( T_e, \; C_t\sqrt{\frac{\nu}{\epsilon}} \right),
\end{equation}

where:

- $\nu = \mu/\rho$ is the molecular kinematic viscosity,
- $C_t$ is a model constant,
- The second term becomes important near walls or when $k/\epsilon$ becomes too large.

If `[scale_limiter=none]`, then $T=T_e$.

### Effective coefficient $C_\mu^{\mathrm{eff}}$

The meaning of $C_\mu$ varies based on the selected turbulence variant:

#### Standard k–$\epsilon$

\begin{equation}
C_\mu^{\mathrm{eff}} = C_\mu.
\end{equation}

#### StandardLowRe

Low-Re models use a damping function $f_\mu$:

\begin{equation}
C_\mu^{\mathrm{eff}} = C_\mu\, f_\mu,
\end{equation}

where $f_\mu$ is computed as:

\begin{equation}
f_\mu = 1 - \exp\!\Big[-(C_{d0}\sqrt{Re_d} + C_{d1}Re_d + C_{d2}Re_d^2)\Big],
\end{equation}

with:

\begin{equation}
Re_d = \frac{\sqrt{k}\,d}{\nu},
\end{equation}

and $d$ is the wall distance (supplied via a functor).
Defaults $C_{d0} = 0.091$, $C_{d1} = 0.0042$, $C_{d2} = 0.00011$ [!cite](lien1996low).

#### Realizable k–$\epsilon$

The realizable model computes $C_\mu$ from the mean strain and rotation
[!cite](shih1995realizable):

\begin{equation}
C_\mu^{\mathrm{eff}} = \frac{1}{A_0 + A_s \, k U^* / \epsilon},
\end{equation}

where:

- $U^* = \sqrt{S_{ij}S_{ij} + \Omega_{ij}\Omega_{ij}}$, with $S_{ij}$ and $\Omega_{ij}$ the mean
  strain-rate and rotation-rate tensors (no system-rotation term),
- $A_s = \sqrt{6}\cos\phi$, $\phi = \tfrac{1}{3}\arccos(\sqrt{6}\,W)$,
  $W = S_{ij}S_{jk}S_{ki}/\tilde{S}^3$, $\tilde{S} = \sqrt{S_{ij}S_{ij}}$,
- $A_0$ is set by [!param](/AuxKernels/kEpsilonViscosity/A0) (default 4.04).

In an equilibrium logarithmic layer this gives $C_\mu^{\mathrm{eff}} \approx 0.09$. This value
replaces $C_\mu$ in the bulk $\mu_t$ only; the wall functions and the two-layer quantities keep the
constant [!param](/AuxKernels/kEpsilonViscosity/C_mu).

!alert warning title=Set $\sigma_\epsilon = 1.2$ for the realizable variants
$\sigma_\epsilon$ is not a parameter of [kEpsilonTKEDSourceSink.md]. It is the `scaling_coeff` of
the $\epsilon$ turbulent diffusion kernel (the one with `diffusion_coeff = mu_t`), or it is built
into a combined coefficient such as $\mu + \mu_t/\sigma_\epsilon$. Its default is 1 and nothing
checks it. `Realizable` and `RealizableTwoLayer` are calibrated with $\sigma_\epsilon = 1.2$ and
$C_2 = 1.9$ [!cite](shih1995realizable); the Standard variants use $\sigma_\epsilon = 1.3$.
`C2_eps` follows the variant automatically, but $\sigma_\epsilon$ has to be set in the input.

#### Two-layer models (StandardTwoLayer and RealizableTwoLayer)

Two-layer models blend between:

- an outer k–$\epsilon$ viscosity,
- an inner near-wall viscosity based on two-layer relations.

Let:

- $\mu_{t,k\epsilon} = \rho C_\mu^{\mathrm{eff}} k T$,
- $\mu_{2L}$ = near-wall turbulent viscosity from the chosen flavor
  (Wolfstein or Norris–Reynolds), see the [theory](theory/turbulenceModeling.md) section for more details.

Using the wall-distance Reynolds number $Re_d = \sqrt{k} d / \nu$,
the two-layer blending function is:

\begin{equation}
\lambda = \tfrac12\left[1 + \tanh\!\left(\frac{Re_d - Re^*}{A}\right)\right],
\end{equation}

with:

- $Re^* = 60$,
- $A$ chosen such that $\lambda \approx 0.98$ at $Re^* + 10$.

The final viscosity is:

\begin{equation}
\mu_t = \lambda\, \mu_{t,k\epsilon} + (1 - \lambda)\, \mu_{2L}.
\end{equation}


## Near-wall bulk treatment

When `[bulk_wall_treatment=true]`, the kernel applies wall functions to compute $\mu_t$ in
wall-bounded cells.

The procedure:

1. Identify the minimum wall distance in the cell.
2. Compute the tangential velocity and friction velocity $u_\tau$.
3. Compute the nondimensional wall distance:

\begin{equation}
y^+ = \frac{\rho u_\tau y}{\mu}.
\end{equation}

4. Classify into:

- *Viscous sublayer*: $y^+ \le 5$ → $\mu_t = 0$,
- *Log‑layer*: $y^+ \ge 30$ → $\mu_t = \mu_\mathrm{log}(y^+)$ (see [theory](theory/turbulenceModeling.md) for more details),
- *Buffer layer*: linear blending between the two regimes.

Four wall-function types are supported:

- `eq_newton`
- `eq_incremental`
- `eq_linearized`
- `neq` (non-equilibrium using local $k$)

Each provides a different method for computing $u_\tau$ and the corresponding $\mu_{wall}$.

This wall-function implementation is fully consistent with the finite-volume INSFVTurbulentViscosity
models used in the Navier–Stokes FV physics.


## Wall distance requirements

The following variants *require* a wall-distance functor:

- `StandardTwoLayer`
- `RealizableTwoLayer`
- `StandardLowRe`

If a required functor is missing, the kernel throws an error during construction.


## Interaction with k and $\epsilon$ kernels

`kEpsilonViscosity` is tightly coupled with the k and $\epsilon$ kernels:

- It uses `k` and $\epsilon$ to compute viscosity.
- The k and $\epsilon$ kernels require $\mu_t$ to compute production and destruction terms.
- Realizable and nonlinear models require strain/rotation invariants also used by k and $\epsilon$ kernels.
- Two-layer models require consistency between:
  - $\epsilon$-equation two-layer region,
  - viscosity two-layer region.

This object must therefore be included in every k–$\epsilon$ turbulence simulation.

!syntax parameters /AuxKernels/kEpsilonViscosity

!syntax inputs /AuxKernels/kEpsilonViscosity

!syntax children /AuxKernels/kEpsilonViscosity
