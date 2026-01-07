# kEpsilonTKEDSourceSink

`kEpsilonTKEDSourceSink` is a finite volume elemental kernel that computes the **turbulent source
and sink terms for the $\epsilon$-equation** in the k–$\epsilon$ family of models. Together with
[`kEpsilonTKESourceSink`](kEpsilonTKESourceSink.md) and [`kEpsilonViscosity`](kEpsilonViscosity.md),
it forms the full set of transport and closure relations for the k–$\epsilon$ turbulence models implemented
in OpenPronghorn.

!alert note
The explanations in this kernel documentation are straightforward.
The reader is referred to the [theory](theory/turbulenceModeling.md) for more details if needed.

The $\epsilon$-equation is written in conservative form as

\begin{equation}
\frac{\partial \rho \epsilon}{\partial t} + \nabla \cdot (\rho \mathbf{u} \epsilon)
= \rho (P_\epsilon - D_\epsilon + S_\epsilon),
\end{equation}

where

- $\epsilon$ is the turbulent dissipation rate,
- $\mathbf{u}$ is the mean velocity,
- $P_\epsilon$ is the production of $\epsilon$,
- $D_\epsilon$ is the destruction term proportional to $\epsilon^2 / k$,
- $S_\epsilon$ collects additional model-specific sources (Yap, low-Re extra production, etc.).

In the implementation, `kEpsilonTKEDSourceSink` forms the $\epsilon$-source term as
\begin{equation}
\text{source}_\epsilon = P_\epsilon - \rho C_{\epsilon2} f_2 \frac{\epsilon^2}{k},
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
production, buoyancy contributions, non-linear production, Yap correction, and low-Re extra
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
- the shear-based $\epsilon$-production term $S_k = \mu_t S^2$ for realizable variants,
- the curvature correction factor $f_c$ (for realizable models when enabled),
- the normalized strain/rotation quantities used by the non-linear constitutive relations.


## $\epsilon$-production term $P_\epsilon$

The $\epsilon$-production term $P_\epsilon$ is built from several contributions:

- $G_k$: shear production of k,
- $G_{\text{nl}}$: non-linear production (optional quadratic/cubic constitutive relations),
- $G'$: low-Re extra production (StandardLowRe only),
- $G_b$: buoyancy production,
- $Y_y$: Yap correction (two-layer and low-Re variants),
- $f_c$: curvature correction factor (realizable variants).

The exact combination depends on the chosen k–$\epsilon$ variant.

### Standard k–$\epsilon$

For the **Standard** high-Re k–$\epsilon$ model, the $\epsilon$-production is

\begin{equation}
P_\epsilon = G_k^{\text{lim}} + G_{\text{nl}} + C_{\epsilon3} G_b,
\end{equation}

where

- $G_k^{\text{lim}}$ is the limited shear production
  \begin{equation}
  G_k^{\text{lim}} = \min(G_k, C_{PL}\,\rho\,\epsilon),
  \end{equation}
- $G_{\text{nl}}$ is the non-linear production term (if `nonlinear_model != none`),
- $C_{\epsilon3}$ is the buoyancy coefficient (usually model-dependent).

The production limiter is consistent with the one used in the k-equation kernel.

### StandardTwoLayer k–$\epsilon$

For the **StandardTwoLayer** model, the $\epsilon$-production is

\begin{equation}
P_\epsilon = G_k + G_{\text{nl}} + C_{\epsilon3} G_b + \rho C_{\epsilon1} Y_y,
\end{equation}

where

- the shear production $G_k$ is taken without limiter,
- the Yap correction $Y_y$ provides an additional near-wall sink/source and is multiplied by
  $\rho C_{\epsilon1}$.

The Yap term is especially important near walls and is detailed in a separate section below.


### StandardLowRe k–$\epsilon$ (low-Re model)

For the **StandardLowRe** model, low-Re corrections are included via damping functions and the
extra production term $G'$. The $\epsilon$-production is

\begin{equation}
P_\epsilon = G_k + G_{\text{nl}} + G' + C_{\epsilon3} G_b + \rho C_{\epsilon1} Y_y,
\end{equation}

where

- $G'$ is the low-Re extra production (see below),
- $f_2$ damping appears in the destruction term,
- Yap correction $Y_y$ is applied similarly to the two-layer model.

### Realizable k–$\epsilon$

For the **Realizable** high-Re model, the $\epsilon$-production is written in terms of the shear-based

production $S_k$:
\begin{equation}
P_\epsilon = f_c S_k + C_{\epsilon3} G_b,
\end{equation}

where

- $S_k = \mu_t S^2$ is the shear contribution,
- $f_c$ is the curvature correction factor (if `curvature_model != none`),
- $C_{\epsilon3} G_b$ is the buoyancy contribution.

In realizable models the non-linear constitutive relation affects the viscosity (via realizable
$C_\mu$) and the k-equation, but $P_\epsilon$ itself is expressed in terms of $S_k$ rather
than $G_k$ or $G_{\text{nl}}$.

### RealizableTwoLayer k–$\epsilon$

The **RealizableTwoLayer** model combines the realizable bulk behavior with a two-layer near-wall
enhancement via Yap:

\begin{equation}
P_\epsilon = f_c S_k + C_{\epsilon3} G_b + \rho C_{\epsilon1} Y_y,
\end{equation}

with the same definitions of $S_k$, $f_c$, and $G_b$ as in the realizable model and the Yap
term as described below.


## Non-linear $\epsilon$-production $G_{\text{nl}}$

When [!param](/LinearFVKernels/kEpsilonTKEDSourceSink/nonlinear_model) is set to `quadratic` or
`cubic`, `kEpsilonTKEDSourceSink` includes the contribution $G_{\text{nl}}$ in the $\epsilon$-production
for Standard-family variants:

- `Standard`:
  \begin{equation}
  P_\epsilon = G_k^{\text{lim}} + G_{\text{nl}} + C_{\epsilon3} G_b,
  \end{equation}

- `StandardTwoLayer`:
  \begin{equation}
  P_\epsilon = G_k + G_{\text{nl}} + C_{\epsilon3} G_b + \rho C_{\epsilon1} Y_y,
  \end{equation}

- `StandardLowRe`:
  \begin{equation}
  P_\epsilon = G_k + G_{\text{nl}} + G' + C_{\epsilon3} G_b + \rho C_{\epsilon1} Y_y.
  \end{equation}

The term $G_{\text{nl}}$ is computed in `TurbulenceMethods` based on the quadratic or cubic
non-linear Reynolds stress models and contracted with $\nabla \mathbf{u}$.


## Curvature correction in $\epsilon$-production

For realizable variants, the curvature correction factor $f_c$ obtained from
`NS::computeCurvatureFactor` is applied multiplicatively to $S_k$ in the $\epsilon$-production:

- `Realizable`:
  \begin{equation}
  P_\epsilon = f_c S_k + C_{\epsilon3} G_b,
  \end{equation}

- `RealizableTwoLayer`:
  \begin{equation}
  P_\epsilon = f_c S_k + C_{\epsilon3} G_b + \rho C_{\epsilon1} Y_y.
  \end{equation}

The same curvature model (e.g. `curvature_model = standard`) is used in both the k-equation and
the $\epsilon$-equation kernels for consistency.


## Yap correction

The **Yap correction** is an additional near-wall term that modifies the $\epsilon$-equation source in
two-layer and low-Re variants. It is enabled when
[!param](/LinearFVKernels/kEpsilonTKEDSourceSink/use_yap) is `true` and a wall-distance functor is
available.

### Yap length scales

The Yap correction is built from two length scales:

- a turbulence length scale $l$,
- an $\epsilon$ length scale $l_\epsilon$.

The implementation uses a *limited* turbulence time scale
to obtain $l$. The effective time scale $T_{\text{eff}}$ is computed as

\begin{equation}
T_1 = \frac{k}{\epsilon},
\end{equation}

\begin{equation}
T_2 = \frac{6\nu}{\epsilon},
\end{equation}

\begin{equation}
T_3 = T_1^{1.625} T_2,
\end{equation}

\begin{equation}
T_{\text{eff}} = T_3^{1/(1.625 + 1)},
\end{equation}

where $\nu = \mu / \rho$ is the kinematic viscosity. The Yap turbulence length scale is then

\begin{equation}
l = \sqrt{T_{\text{eff}} \epsilon}.
\end{equation}

The $\epsilon$ length scale is defined as

\begin{equation}
l_\epsilon = C_l d,
\end{equation}

where $d$ is the wall distance and

\begin{equation}
C_l = 0.42 C_\mu^{-3/4}
\end{equation}

with $C_\mu$ the model’s closure constant (or a realizable value when applicable).

### Yap correction term $Y_y$

The Yap correction term $\gamma_y$ is computed in `TurbulenceMethods::computeGammaY` as

\begin{equation}
\gamma_y =
C_w \frac{\epsilon^2}{k}
\max\left[\left(\frac{l}{l_\epsilon} - 1\right)\left(\frac{l}{l_\epsilon}\right)^2, 0\right],
\end{equation}

where $C_w$ is a model coefficient. In the $\epsilon$-equation, this is translated into a source term

\begin{equation}
\rho C_{\epsilon1} Y_y = \rho C_{\epsilon1} \gamma_y,
\end{equation}

and added to $P_\epsilon$ for the variants that include Yap.

The Yap contribution appears in:

- `StandardTwoLayer`,
- `StandardLowRe`,
- `RealizableTwoLayer`.


## Low-Re extra production $G'$

For the **StandardLowRe** model, an additional low-Re production term $G'$ is included in
$P_\epsilon$ when [!param](/LinearFVKernels/kEpsilonTKEDSourceSink/use_low_re_Gprime) is `true`
and a wall-distance functor is provided.

The term $G'$ is defined as

\begin{equation}
G' = D f_2 \left(G_k + 2 \mu_t \frac{k}{d^2}\right) \exp(-E \mathrm{Re}_d^2),
\end{equation}

where

- $D$ and $E$ are model coefficients,
- $f_2$ is the low-Re damping function,
  \begin{equation}
  f_2 = 1 - C \exp(-\mathrm{Re}_t^2),
  \end{equation}
  with $C$ a low-Re coefficient,
- $\mathrm{Re}_d = \sqrt{k} d / \nu$ is the wall-distance Reynolds number,
- $\mathrm{Re}_t = k^2 / (\nu \epsilon)$ is the turbulence Reynolds number.

In the $\epsilon$-production for the StandardLowRe variant,

\begin{equation}
P_\epsilon = G_k + G_{\text{nl}} + G' + C_{\epsilon3} G_b + \rho C_{\epsilon1} Y_y.
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

- `Standard`, `StandardTwoLayer`, `Realizable`, `RealizableTwoLayer`:
  \begin{equation}
  f_2 = 1,
  \end{equation}

- `StandardLowRe`:
  \begin{equation}
  f_2 = 1 - C \exp(-\mathrm{Re}_t^2),
  \end{equation}
  consistent with the low-Re formulation used for $G'$.

The same `f2_SKE_LRe` function is used for both the extra production term and the destruction
damping in the low-Re model.


!syntax parameters /LinearFVKernels/kEpsilonTKEDSourceSink

!syntax inputs /LinearFVKernels/kEpsilonTKEDSourceSink

!syntax children /LinearFVKernels/kEpsilonTKEDSourceSink
