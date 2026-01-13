# kEpsilonTKESourceSink

`kEpsilonTKESourceSink` is a finite volume elemental kernel that computes the *turbulent source
and sink term for the turbulent kinetic energy equation* in the k–$\epsilon$ family of models.

!alert note
The explanations in this kernel documentation are straightforward.
The reader is referred to the [theory](theory/turbulenceModeling.md) for more details if needed.

The turbulent kinetic energy (TKE) equation is written in conservative form as

\begin{equation}
\frac{\partial \rho k}{\partial t} + \nabla \cdot (\rho \mathbf{u} k)
= \rho (P_k + G_b + G_{\text{nl}} - \gamma_M - \epsilon),
\end{equation}

where:

- $k$ is the turbulent kinetic energy,
- $\mathbf{u}$ is the mean velocity,
- $P_k$ is the shear production,
- $G_b$ is the buoyancy production,
- $G_{\text{nl}}$ is the production due to non-linear Reynolds stresses,
- $\gamma_M$ is a compressibility correction, and
- $\epsilon$ is the dissipation rate supplied by the $\epsilon$-equation.

`kEpsilonTKESourceSink` forms the *right-hand side source term* for the k-equation by computing
the combination

\begin{equation}
\text{source} = P_k + G_b + G_{\text{nl}} - \gamma_M - \epsilon,
\end{equation}

with expressions that depend on the chosen k–$\epsilon$ variant and options.

The kernel is designed to be used together with:

- [`kEpsilonTKEDSourceSink`](kEpsilonTKEDSourceSink.md) for the $\epsilon$-equation, and
- [`kEpsilonViscosity`](kEpsilonViscosity.md) for the turbulent viscosity $\mu_t$.

## Model variants

The kernel supports several members of the k–$\epsilon$ family via the
[!param](/LinearFVKernels/kEpsilonTKESourceSink/k_epsilon_variant) parameter:

- `Standard` — classical high-Reynolds-number k–$\epsilon$ model.
- `StandardLowRe` — low-Reynolds-number k–$\epsilon$ model with damping functions.
- `StandardTwoLayer` — two-layer k–$\epsilon$ formulation using a blending between outer k–$\epsilon$ and
  near-wall length scales.
- `Realizable` — realizable k–$\epsilon$ model with variable $C_\mu$ depending on strain/rotation
  invariants.
- `RealizableTwoLayer` — realizable model with a two-layer near-wall treatment.

The choice of variant affects:

- how the shear production is formed,
- whether a curvature correction factor is applied,
- how near-wall behavior is modeled (two-layer vs. high-Re),
- whether low-Re damping functions and extra terms are used (via $\epsilon$-equation).


## Bulk formulation

In *bulk cells* (not explicitly treated as near-wall by this kernel), the source term is
computed in `computeBulkProduction()` as

\begin{equation}
\text{source}_{\text{bulk}} = P_k + G_b + G_{\text{nl}} - \gamma_M - \epsilon.
\end{equation}

The individual contributions are detailed below.

### Strain-rate tensor and invariants

The mean velocity gradient $\nabla \mathbf{u}$ is used to form the symmetric strain-rate tensor
$\mathbf{S}$ and the antisymmetric rotation tensor $\mathbf{W}$:

\begin{equation}
\mathbf{S} = \frac{1}{2} \left(\nabla \mathbf{u} + \nabla \mathbf{u}^T \right)
\end{equation}

\begin{equation}
\mathbf{W} = \frac{1}{2} \left(\nabla \mathbf{u} - \nabla \mathbf{u}^T \right)
\end{equation}

`kEpsilonTKESourceSink` calls
`NS::computeStrainRotationInvariants(u, v, w, elem_arg, state)` from
`TurbulenceMethods` to obtain the invariants:

- $S^2 = 2 S_{ij} S_{ij}$,
- $W^2 = 2 W_{ij} W_{ij}$,
- $\nabla \cdot \mathbf{u}$.

These invariants are shared by the turbulence viscosity, non-linear models, and curvature
correction.

### Shear production $G_k$

The base formulation of the shear production is

\begin{equation}
G_k = \mu_t S^2,
\end{equation}

where

- $\mu_t$ is the turbulent viscosity provided as a functor, and
- $S = \sqrt{2 \mathbf{S}:\mathbf{S}} = \sqrt{S^2}$.

For compressible flows, the implementation allows inclusion of the extra trace terms

\begin{equation}
G_k = \mu_t S^2
- \frac{2}{3}\left(\rho k \nabla \cdot \mathbf{u}
                 + \mu_t (\nabla\cdot\mathbf{u})^2\right),
\end{equation}

consistent with the constitutive relation for the Reynolds stress tensor in compressible flow.

The shear production is then modified depending on the k–$\epsilon$ variant and curvature correction model
(see below).

### Curvature / rotation correction $f_c$

For realizable variants, an optional curvature correction factor $f_c$ may be applied via
[!param](/LinearFVKernels/kEpsilonTKESourceSink/curvature_model). The current implementation
supports:

- `curvature_model = none` — no curvature correction, $f_c = 1$;
- `curvature_model = standard` — Spalart–Shur–type rotation/curvature correction based on
  $S^2$ and $W^2$.

`NS::computeCurvatureFactor` builds a correction $f_c = f_c(S^2, W^2)$ that enhances production
in stabilizing curvature and reduces it in destabilizing regions. The corrected shear production is

- *Standard*, *StandardLowRe*, *StandardTwoLayer*:
  \begin{equation}
  P_k = G_k;
  \end{equation}

- *Realizable*, *RealizableTwoLayer*:
  \begin{equation}
  P_k = f_c G_k.
  \end{equation}

### Buoyancy production $G_b$

If [!param](/LinearFVKernels/kEpsilonTKESourceSink/use_buoyancy) is `true`, buoyancy production is
added via

\begin{equation}
G_b = \beta \frac{\mu_t}{Pr_t} \left(\nabla T \cdot \mathbf{g}\right),
\end{equation}

where

- $\beta$ is the thermal expansion coefficient,
- $Pr_t$ is the turbulent Prandtl number,
- $T$ is the temperature,
- $\mathbf{g}$ is the gravitational acceleration.

The gradient $\nabla T$ and gravity vector are supplied as functors / parameters, and the
computation is performed through `NS::computeGb`. Buoyancy can either augment or reduce total
production depending on the sign of the dot product $\nabla T\cdot\mathbf{g}$.

### Non-linear production $G_{\text{nl}}$

`kEpsilonTKESourceSink` supports optional *non-linear constitutive relations* for the Reynolds
stress tensor. The choice is controlled by
[!param](/LinearFVKernels/kEpsilonTKESourceSink/nonlinear_model):

- `nonlinear_model = none` — no non-linear production, $G_{\text{nl}} = 0$;
- `nonlinear_model = quadratic` — quadratic constitutive relation;
- `nonlinear_model = cubic` — quadratic + cubic terms.

The non-linear Reynolds stress contribution is computed in `TurbulenceMethods` as

\begin{equation}
G_{\text{nl}} = \mathbf{T}_{\text{RANS,NL}} : \nabla \mathbf{u},
\end{equation}

where $\mathbf{T}_{\text{RANS,NL}}$ is obtained from `NS::computeTRANS_NL` and depends on:

- the velocity gradient $\nabla \mathbf{u}$ (through $\mathbf{S}$ and $\mathbf{W}$),
- the invariants $S^2$, $W^2$,
- $\mu_t$, $k$, $\epsilon$,
- fixed model coefficients $C_{\text{NL}1\ldots 7}$ and a variable $C_\mu$ computed by
  `NS::Cmu_nonlinear`.

In all k–$\epsilon$ variants, `kEpsilonTKESourceSink` adds $G_{\text{nl}}$ to the k-equation’s production:

\begin{equation}
\text{source} = (P_k + G_b + G_{\text{nl}} - \gamma_M) - \epsilon.
\end{equation}

### Compressibility correction $\gamma_M$

If [!param](/LinearFVKernels/kEpsilonTKESourceSink/use_compressibility) is `true`, a compressibility
correction is added to the sink term:

\begin{equation}
\gamma_M = \rho C_M \frac{k \epsilon}{c^2},
\end{equation}

where $C_M$ is a model coefficient and $c$ is the speed of sound. This term is computed by
`NS::computeGammaM` and is subtracted from total production, representing the transfer of TKE into
acoustic modes in compressible flow.

### Production limiter

To avoid excessive TKE generation in stagnation regions, a *production limiter* is applied:

\begin{equation}
G_k \le C_{PL}\,\rho\,\epsilon.
\end{equation}

The limiter coefficient $C_{PL}$ is available as
[!param](/LinearFVKernels/kEpsilonTKESourceSink/C_pl). This limiter follows the standard approach
from Durbin (1996) and Menter (1994) and is applied in the bulk region.

### Summary of bulk source

In non-wall-bounded cells, the kernel assembles the following source term for the residual:

\begin{equation}
\text{source}_{\text{bulk}}
= \underbrace{P_k}_{\text{shear (varies by variant)}}
+ \underbrace{G_b}_{\text{buoyancy}}
+ \underbrace{G_{\text{nl}}}_{\text{non-linear model}}
- \underbrace{\gamma_M}_{\text{compressibility}}
- \underbrace{\epsilon}_{\text{dissipation}}.
\end{equation}


## Wall formulation

Cells adjacent to boundaries listed in
[!param](/LinearFVKernels/kEpsilonTKESourceSink/walls) are treated as *next-to-wall cells* with a
different expression for production and destruction, similar in spirit to
[`LinearFVTKESourceSink`](LinearFVTKESourceSink.md). Wall distances are typically computed using
[`WallDistanceAux`](WallDistanceAux.md), which returns the distance from each cell centroid to the
nearest wall.

The near-wall treatment is implemented in the matrix contribution for this kernel, where a
distinction is made between the *viscous sublayer* and the *logarithmic region* based on the
non-dimensional wall distance

\begin{equation}
y^+ = \frac{\rho y_p u_\tau}{\mu},
\end{equation}

with

- $y_p$ the distance from the cell center to the wall,
- $u_\tau$ the friction velocity obtained from a wall-function model,
- $\mu$ the molecular viscosity.

The same wall-function options as used for the turbulent viscosity are available via
[!param](/LinearFVKernels/kEpsilonTKESourceSink/wall_treatment) (e.g. `eq_newton`,
`eq_incremental`, `eq_linearized`, `neq`).

### Viscous sublayer ($y^+ < 11.25$)

For cells in the viscous sublayer (sub-laminar region), turbulent production is negligible and the
model sets

\begin{equation}
G_k = 0.
\end{equation}

Destruction is modeled using a near-wall expression of the form

\begin{equation}
\epsilon = \frac{2\mu k}{y_p^2},
\end{equation}
which enforces rapid dissipation very close to the wall.

### Logarithmic region ($y^+ \ge 11.25$)

For cells in the log layer, production and destruction take the classical wall-function forms:

- Production:
  \begin{equation}
  G_k = \tau_w \|\nabla \mathbf{u}\|
      = \mu_w \|\nabla \mathbf{u}\| \frac{C_\mu^{1/4}\sqrt{k}}{\kappa y_p},
  \end{equation}
  where
  - $\tau_w$ is the wall shear stress,
  - $\mu_w$ is the effective wall viscosity (including turbulent contributions),
  - $\|\nabla \mathbf{u}\|$ is the wall-normal velocity gradient norm,
  - $C_\mu$ is the k–$\epsilon$ closure coefficient,
  - $\kappa \approx 0.41$ is the von Kármán constant.

- Destruction:
  \begin{equation}
  \epsilon = C_\mu^{3/4} \frac{\rho k^{3/2}}{\kappa y_p}.
  \end{equation}

In the logarithmic region, production and destruction are balanced in such a way that the near-wall
TKE profile is consistent with the log-law of the wall.

!alert note
When the near-wall treatment is handled by this kernel, the wall boundary condition for $k$
is effectively provided by the model and explicit Dirichlet boundary conditions for $k$ on
those walls are typically unnecessary.

## Interaction with $\epsilon$-equation and viscosity

`kEpsilonTKESourceSink` only forms the *k-equation* source terms. For a complete k–$\epsilon$ model, this
kernel must be used together with:

- [`kEpsilonTKEDSourceSink`](kEpsilonTKEDSourceSink.md), which provides the $\epsilon$-equation sources,
  including:
  - realizable production terms,
  - non-linear contributions,
  - low-Re extra production $G'$,
  - Yap correction,
  - buoyancy source for $\epsilon$;
- [`kEpsilonViscosity`](kEpsilonViscosity.md), which computes the turbulent viscosity $\mu_t$
  from the current $k$ and $\epsilon$, including:
  - low-Re damping for SKE-LRe,
  - two-layer viscosity blending,
  - realizable variable $C_\mu$,
  - bulk wall treatment for $\mu_t$ (if requested).

The three objects together implement the full k–$\epsilon$ turbulence model family.`

!syntax parameters /LinearFVKernels/kEpsilonTKESourceSink

!syntax inputs /LinearFVKernels/kEpsilonTKESourceSink

!syntax children /LinearFVKernels/kEpsilonTKESourceSink
