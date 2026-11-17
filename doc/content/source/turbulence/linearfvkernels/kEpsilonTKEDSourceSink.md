# kEpsilonTKEDSourceSink

`kEpsilonTKEDSourceSink` is a finite volume elemental kernel that computes the *turbulent source
and sink terms for the $\epsilon$-equation* in the k–$\epsilon$ family of models. Together with
[`kEpsilonTKESourceSink`](kEpsilonTKESourceSink.md) and [`kEpsilonViscosity`](kEpsilonViscosity.md),
it forms the full set of transport and closure relations for the k–$\epsilon$ turbulence models implemented
in OpenPronghorn.

!alert note
The explanations in this kernel documentation are a summary.
The reader is referred to the [theory](theory/turbulenceModeling.md) for more details if needed.

The $\epsilon$-equation is written in conservative form as

\begin{equation}
\frac{\partial \rho \epsilon}{\partial t} + \nabla \cdot (\rho \mathbf{u} \epsilon)
= P_\epsilon - D_\epsilon,
\end{equation}

where

- $\epsilon$ is the turbulent dissipation rate,
- $\mathbf{u}$ is the mean velocity,
- $P_\epsilon$ is the *assembled production-side term* for the $\epsilon$-equation (shear, buoyancy,
  nonlinear production, and any variant-specific corrections such as Yap or the low-Re $G'$ term),
- $D_\epsilon$ is the destruction term proportional to $\rho\,\epsilon^2/k$.

Some references split the right-hand side as $\rho\,(P_\epsilon - D_\epsilon + S_\epsilon)$. In
OpenPronghorn we fold any additional sources/corrections that would fall under $S_\epsilon$ directly into
$P_\epsilon$, so the kernel assembles a single production term.

In the implementation, `kEpsilonTKEDSourceSink` forms the $\epsilon$-source term as
\begin{equation}
\text{source}_\epsilon = P_\epsilon - D_\epsilon
= P_\epsilon - \rho C_{\epsilon2} f_2 \frac{\epsilon^2}{k},
\end{equation}
where $C_{\epsilon2}$ and $f_2$ depend on the selected k–$\epsilon$ variant and damping model.


## Model variants

The kernel supports the same k–$\epsilon$ variants as [kEpsilonTKESourceSink.md] via the
[!param](/LinearFVKernels/kEpsilonTKEDSourceSink/k_epsilon_variant) parameter:

- `Standard` — classical high-Reynolds-number k–$\epsilon$ model.
- `StandardLowRe` — low-Reynolds-number k–$\epsilon$ model with damping functions.
- `StandardTwoLayer` — two-layer k–$\epsilon$ formulation using a blending between outer k–$\epsilon$ and
  near-wall length scales.
- `Realizable` — realizable k–$\epsilon$ model with variable $C_\mu$.
- `RealizableTwoLayer` — realizable model with a two-layer near-wall treatment.

The $\epsilon$-production term $P_\epsilon$ is assembled differently for each variant, combining shear
production, buoyancy contributions, nonlinear production, Yap correction, and low-Re extra
production $G'$.


## Strain-rate invariants and basic quantities

As in the k-equation kernel, `kEpsilonTKEDSourceSink` uses the invariants of the strain and rotation
tensors:

- $S^2 = 2 S_{ij} S_{ij}$,
- $W^2 = 2 W_{ij} W_{ij}$,
- $\nabla \cdot \mathbf{u}$,

See the [theory](theory/turbulenceModeling.md) section for the definition of these ones.
These invariants are used to construct:

- the shear production of k,
  \begin{equation}
  G_k = \mu_t S^2 \;\; (\text{plus compressible terms if enabled}),
  \end{equation}
  With [!param](/LinearFVKernels/kEpsilonTKEDSourceSink/use_compressibility) set to `true`, the
  $\tfrac{2}{3}$ dilatation terms are included in $G_k$ (set the same value as in the $k$ kernel).
  The compressibility sink $\gamma_M$ acts in the $k$ equation only (see
  [kEpsilonTKESourceSink.md]), so this kernel needs no `speed_of_sound` or `C_M`.
- the strain magnitude $S = \sqrt{S^2}$ and $\eta = S k / \epsilon$ used by the realizable
  production $\rho C_1 S \epsilon$,
- the curvature correction factor $f_c$ (for realizable models when enabled),
- the normalized strain/rotation quantities used by the nonlinear constitutive relations.


## $\epsilon$-production term $P_\epsilon$

The $\epsilon$-production term $P_\epsilon$ is built from several contributions. Here, $P_\epsilon$ denotes the full non-destruction part of the $\epsilon$ source term (i.e., it already includes terms that are sometimes written separately as $S_\epsilon$, such as low-Re corrections). The Yap source $\rho\gamma_y$ is added separately (see below):

- $G_k$: shear production of k,
- $G_{\text{nl}}$: nonlinear production (optional quadratic/cubic constitutive relations),
- $G'$: low-Re extra production (StandardLowRe only),
- $G_b$: buoyancy production,
- $f_c$: curvature correction factor (realizable variants).

The exact combination depends on the chosen k–$\epsilon$ variant.

### Standard k–$\epsilon$

For the *Standard* high-Re k–$\epsilon$ model, the $\epsilon$-production is

\begin{equation}
G_k^{\text{lim}} + G_{\text{nl}} + C_{\epsilon3} G_b,
\end{equation}

where

- $G_k^{\text{lim}}$ is the limited shear production
  \begin{equation}
  G_k^{\text{lim}} = \min(G_k, C_{PL}\,\rho\,\epsilon),
  \end{equation}
- $G_{\text{nl}}$ is the nonlinear production term (if `nonlinear_model != none`),
- $C_{\epsilon3}$ is the buoyancy coefficient (usually model-dependent).

When `use_buoyancy = true` and the gravity vector is nonzero, $C_{\epsilon3}$ is computed in each
cell as

\begin{equation}
C_{\epsilon3} = \tanh\left|\frac{u_\parallel}{u_\perp}\right|,
\end{equation}

where $u_\parallel$ and $u_\perp$ are the velocity components parallel and perpendicular to
$\mathbf{g}$. The coefficient is close to 1 in vertical buoyant shear layers and close to 0 in
horizontal ones. The user value `C3_eps` is used only where $u_\perp = 0$. The same $C_{\epsilon3} G_b$
term is used by every variant, including the realizable ones, where it enters as
$C_{\epsilon1} C_{\epsilon3} G_b\,\epsilon/k$.

The production limiter is consistent with the one used in the k-equation kernel.

### StandardTwoLayer k–$\epsilon$

For the *StandardTwoLayer* model, the $\epsilon$-production is

\begin{equation}
G_k + G_{\text{nl}} + C_{\epsilon3} G_b,
\end{equation}

where the shear production $G_k$ is taken without limiter. The optional Yap source is described
below.


### StandardLowRe k–$\epsilon$ (low-Re model)

For the *StandardLowRe* model, low-Re corrections are included via damping functions and the
extra production term $G'$. The $\epsilon$-production is

\begin{equation}
G_k + G_{\text{nl}} + G' + C_{\epsilon3} G_b,
\end{equation}

where

- $G'$ is the low-Re extra production (see below),
- $f_2$ damping appears in the destruction term,
- the optional Yap source is the same as in the two-layer model.

### Realizable k–$\epsilon$

For the *Realizable* high-Re model [!cite](shih1995realizable), the shear production of
$\epsilon$ does not contain $G_k$. The $\epsilon$ source is

\begin{equation}
f_c \, \rho C_1 S \epsilon + C_{\epsilon1} \frac{\epsilon}{k} C_{\epsilon3} G_b,
\qquad
C_1 = \max\left(0.43, \frac{\eta}{\eta + 5}\right), \qquad \eta = S \frac{k}{\epsilon},
\end{equation}

where

- $S = \sqrt{2 S_{ij} S_{ij}}$ is the strain magnitude,
- $f_c$ is the curvature correction factor (if `curvature_model != none`),
- the shear term $\rho C_1 S \epsilon$ is added directly, not through $C_{\epsilon1}/T_e$,
- only the buoyancy term uses the fixed $C_{\epsilon1}$ = [!param](/LinearFVKernels/kEpsilonTKEDSourceSink/C1_eps).

$G_{\text{nl}}$ and $G'$ do not enter the realizable $\epsilon$ equation.

### RealizableTwoLayer k–$\epsilon$

The *RealizableTwoLayer* model uses the realizable bulk $\epsilon$ source above, blended with the
two-layer near-wall $\epsilon$. The Yap correction is not available for this variant.


## nonlinear $\epsilon$-production $G_{\text{nl}}$

When [!param](/LinearFVKernels/kEpsilonTKEDSourceSink/nonlinear_model) is set to `quadratic` or
`cubic`, `kEpsilonTKEDSourceSink` includes the contribution $G_{\text{nl}}$ in the $\epsilon$-production
for the Standard-family variants (selecting a nonlinear relation with a realizable variant is an error):

- `Standard`:
  \begin{equation}
  \text{Production}= G_k^{\text{lim}} + G_{\text{nl}} + C_{\epsilon3} G_b,
  \end{equation}

- `StandardTwoLayer`:
  \begin{equation}
  \text{Production} = G_k + G_{\text{nl}} + C_{\epsilon3} G_b,
  \end{equation}

- `StandardLowRe`:
  \begin{equation}
  \text{Production} = G_k + G_{\text{nl}} + G' + C_{\epsilon3} G_b.
  \end{equation}

The term $G_{\text{nl}}$ is based on the quadratic or cubic
nonlinear Reynolds stress models and contracted with $\nabla \mathbf{u}$.


## Curvature correction in $\epsilon$-production

For realizable variants, the curvature correction factor $f_c$
is applied multiplicatively to $S_k$ in the $\epsilon$-production:

- `Realizable`:
  \begin{equation}
  \text{Production} = f_c S_k + C_{\epsilon3} G_b,
  \end{equation}

- `RealizableTwoLayer`:
  \begin{equation}
  \text{Production} = f_c S_k + C_{\epsilon3} G_b.
  \end{equation}

The same curvature model (e.g. `curvature_model = standard`) is used in both the k-equation and
the $\epsilon$-equation kernels for consistency.


## Yap correction

The Yap correction [!cite](yap1987correction) is enabled with
[!param](/LinearFVKernels/kEpsilonTKEDSourceSink/use_yap) and is only available for the
`StandardTwoLayer` and `StandardLowRe` variants (any other variant is an error). It needs the
wall distance $d$. The term adds the source $\rho\gamma_y$ directly to the $\epsilon$ equation
(it does not pass through $C_{\epsilon1}/T_e$):

\begin{equation}
\gamma_y =
C_w \frac{\epsilon^2}{k}
\max\left[\left(\frac{l}{l_\epsilon} - 1\right)\left(\frac{l}{l_\epsilon}\right)^2, 0\right],
\qquad l = \frac{k^{3/2}}{\epsilon}, \qquad l_\epsilon = C_l d,
\end{equation}

with $C_w = 0.83$ ([!param](/LinearFVKernels/kEpsilonTKEDSourceSink/Cw)) and
$C_l = 0.42 C_\mu^{-3/4} \approx 2.55$ (from the constant
[!param](/LinearFVKernels/kEpsilonTKEDSourceSink/C_mu)). In an equilibrium log layer $l = l_\epsilon$,
so the term is zero. It is non-negative and acts where the turbulence length exceeds its
equilibrium value, e.g. near reattachment or impingement. In two-layer variants it is part of the
bulk source and is blended with $\lambda$.


## Low-Re extra production $G'$

For the *StandardLowRe* model, an additional low-Re production term $G'$ is included in
$P_\epsilon$ when [!param](/LinearFVKernels/kEpsilonTKEDSourceSink/use_low_re_Gprime) is `true`
and a wall-distance functor is provided.

The term $G'$ is defined as

\begin{equation}
G' = D f_2 \left(G_k + 2 \mu \frac{k}{d^2}\right) \exp(-E \mathrm{Re}_d^2),
\end{equation}

where

- $D$ and $E$ are model coefficients,
- $f_2$ is the low-Re damping function,
  \begin{equation}
  f_2 = 1 - C \exp(-\mathrm{Re}_t^2),
  \end{equation}
  with $C$ a low-Re coefficient,
- $\mu$ is the molecular dynamic viscosity, so that $2\mu k/d^2 = \rho\epsilon_w$ is the wall
  value of $\epsilon$,
- $\mathrm{Re}_d = \sqrt{k} d / \nu$ is the wall-distance Reynolds number,
- $\mathrm{Re}_t = k^2 / (\nu \epsilon)$ is the turbulence Reynolds number.

Defaults: $C = 0.3$ (`C_lowRe`), $D = 1$ (`D_lowRe`), $E = 0.00375$ (`E_lowRe`)
[!cite](lien1996low).

In the $\epsilon$-production for the StandardLowRe variant the production term is

\begin{equation}
G_k + G_{\text{nl}} + G' + C_{\epsilon3} G_b.
\end{equation}

The term $G'$ vanishes as $\mathrm{Re}_d \to \infty$ due to the exponential damping and is
most active in the near-wall low-Re region.


## Destruction term and damping

The destruction term in the $\epsilon$-equation is modeled as

\begin{equation}
D_\epsilon = \rho C_{\epsilon2} f_2 \frac{\epsilon^2}{k},
\end{equation}

where $C_{\epsilon2}$ is the model coefficient and $f_2$ is a damping function that depends on
the k–$\epsilon$ variant:

- `Standard`, `StandardTwoLayer`:
  \begin{equation}
  f_2 = 1,
  \end{equation}

- `Realizable`, `RealizableTwoLayer` [!cite](shih1995realizable):
  \begin{equation}
  f_2 = \frac{k}{k + \sqrt{\nu\epsilon}},
  \end{equation}
  so that $D_\epsilon = \rho C_{\epsilon2} \epsilon^2 / (k + \sqrt{\nu\epsilon})$, which stays finite
  as $k \to 0$,

- `StandardLowRe`:
  \begin{equation}
  f_2 = 1 - C \exp(-\mathrm{Re}_t^2),
  \end{equation}
  consistent with the low-Re formulation used for $G'$.

The same function is used for both the extra production term and the destruction
damping in the low-Re model.

### Model constants

The default of [!param](/LinearFVKernels/kEpsilonTKEDSourceSink/C2_eps) depends on the variant:
1.92 for the Standard variants [!cite](launder1974numerical) and 1.9 for the Realizable variants
[!cite](shih1995realizable). A value set in the input is always used.

The turbulent Prandtl number $\sigma_\epsilon$ is not a parameter of this kernel: it is the
`scaling_coeff` of the $\epsilon$ turbulent diffusion kernel (`LinearFVTurbulentDiffusion` or
`LinearFVTurbulentBlendedDiffusion`), whose default is 1. When
[!param](/LinearFVKernels/kEpsilonTKEDSourceSink/check_sigma_eps) is `true` (default), the kernel
warns if a constant `scaling_coeff` of an $\epsilon$ diffusion kernel differs from the calibrated
value of the selected variant: $\sigma_\epsilon = 1.3$ for the Standard variants and
$\sigma_\epsilon = 1.2$ for the Realizable variants. A diffusion kernel without `scaling_coeff`
(for example the molecular `LinearFVDiffusion`) and a non-constant functor are not checked.


!syntax parameters /LinearFVKernels/kEpsilonTKEDSourceSink

!syntax inputs /LinearFVKernels/kEpsilonTKEDSourceSink

!syntax children /LinearFVKernels/kEpsilonTKEDSourceSink
