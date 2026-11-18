# $k$–$\epsilon$ Turbulence Modeling in OpenPronghorn

This page documents the Reynolds-averaged Navier–Stokes (RANS) turbulence-modeling
capabilities implemented in OpenPronghorn. It focuses on the $k$–$\epsilon$ family of
eddy-viscosity models and practical guidance for choosing and using each model in applications.

For general background on turbulence and RANS closures, see Pope [!cite](pope2000turbulent)
and Wilcox [!cite](wilcox2006turbulence). The baseline $k$–$\epsilon$ model follows
Launder and Spalding [!cite](launder1974numerical), with additional options inspired by
realizable, low-Reynolds-number, two-layer, curvature/rotation, buoyancy, and
compressibility corrections [!cite](shih1995realizable,spalart1997rotation,yap1987correction,wolfstein1969velocity,norris1975oneequation,jongen1998simulation).


## Introduction

Most OpenPronghorn applications use RANS turbulence modeling to represent the effect of
unresolved turbulent fluctuations on the mean flow. The main motivations are:

- *Cost*: Direct numerical simulation is generally infeasible at reactor-relevant Reynolds
  numbers, and even LES can be prohibitively expensive for large, geometrically complex
  configurations.
- *Robustness*: Two-equation eddy-viscosity models are well established, robust, and
  widely used in engineering practice.
- *Integration with thermal–hydraulics*: The use of RANS models fits naturally into
  finite-volume formulations for coupled momentum, energy, and species transport equations.

OpenPronghorn’s turbulence infrastructure is designed to be:

- *Modular*: The same utilities are reused by several $k$–$\epsilon$ kernels and auxiliary objects.
- *Extensible*: New turbulence models can be added by leveraging the existing functor and
  finite-volume infrastructure.
- *Configurable*: Users can select a model variant and enable or disable individual
  corrections (low-Re, two-layer, Yap, nonlinear, curvature, buoyancy, compressibility)
  through input-file parameters.


## Theory of RANS modeling

### Reynolds decomposition and RANS equations

Let $\mathbf{u}(\mathbf{x},t)$ and $p(\mathbf{x},t)$ denote the instantaneous
velocity and pressure of a constant-property Newtonian fluid with density $\rho$ and
dynamic viscosity $\mu$. In RANS modeling we write

\begin{equation}
\mathbf{u} = \overline{\mathbf{u}} + \mathbf{u}', \qquad
p = \overline{p} + p',
\end{equation}

where the overbar denotes a suitable average (e.g., time or ensemble) and primes denote
fluctuations [!cite](pope2000turbulent).

Substituting into the Navier–Stokes equations and averaging yields the mean-momentum
equations

\begin{equation}
\rho \left( \frac{\partial \overline{\mathbf{u}}}{\partial t}
+ \overline{\mathbf{u}} \cdot \nabla \overline{\mathbf{u}} \right)
= - \nabla \overline{p}
+ \nabla \cdot
  \left[ \mu \left( \nabla \overline{\mathbf{u}} +
  (\nabla \overline{\mathbf{u}})^T \right) \right]
- \nabla \cdot \boldsymbol{\tau}^R
+ \rho \mathbf{g},
\end{equation}

with the *Reynolds-stress tensor*

\begin{equation}
\boldsymbol{\tau}^R = - \rho \, \overline{\mathbf{u}' \mathbf{u}'}.
\end{equation}

The appearance of $\boldsymbol{\tau}^R$ introduces the *closure problem* of turbulence:
additional equations or models are required to express the Reynolds stresses in terms of
mean quantities.

### Boussinesq eddy-viscosity hypothesis

In classical eddy-viscosity models we assume that the anisotropic part of the Reynolds
stresses is proportional to the mean strain-rate tensor

\begin{equation}
\boldsymbol{\tau}^R =
2 \mu_t \mathbf{S}
- \frac{2}{3} \rho k \mathbf{I},
\end{equation}

where

- $\mu_t$ is the turbulent (eddy) viscosity,
- $k = \tfrac{1}{2} \overline{u'_i u'_i}$ is the turbulent kinetic energy,
- $\mathbf{S}$ is the strain-rate tensor with components
  \begin{equation}
  S_{ij} = \frac{1}{2}
  \left( \frac{\partial \overline{u_i}}{\partial x_j} +
         \frac{\partial \overline{u_j}}{\partial x_i} \right).
  \end{equation}

In OpenPronghorn, $\mu_t$ is provided by the $k$–$\epsilon$ models via
[`kEpsilonViscosity`](kEpsilonViscosity.md). The Boussinesq hypothesis is reasonably accurate for many wall-bounded
shear flows, but is known to underperform in strongly anisotropic, rotating, or
curvature-dominated flows [!cite](pope2000turbulent,wilcox2006turbulence). To mitigate these
limitations OpenPronghorn includes:

- optional *nonlinear stress* models (quadratic and cubic) that extend the constitutive
  relation beyond linear dependence on $\mathbf{S}$;
- a *curvature/rotation correction* that modifies production based on invariants of
  $\mathbf{S}$ and the rotation tensor $\mathbf{W}$;
- *realizable* formulations in which $C_\mu$ depends on local invariants of
  $\mathbf{S}$ and $\mathbf{W}$ [!cite](shih1995realizable).


## The $k$–$\epsilon$ models implemented

### General transport equations

The baseline OpenPronghorn implementation follows the two-equation model of
Launder and Spalding [!cite](launder1974numerical). In conservative form the transport
equations for $k$ and $\epsilon$ are

\begin{equation}
\frac{\partial \rho k}{\partial t} + \nabla \cdot (\rho \mathbf{u} k)
=  \nabla \cdot \left[ \rho
    \left( \nu + \frac{\nu_t}{\sigma_k} \right) \nabla k
  \right]
 + P_k + G_b - \rho \epsilon + \rho S_k,
\end{equation}

\begin{equation}
\frac{\partial \rho \epsilon}{\partial t} + \nabla \cdot (\rho \mathbf{u} \epsilon)
=  \nabla \cdot \left[ \rho
    \left( \nu + \frac{\nu_t}{\sigma_\epsilon} \right) \nabla \epsilon
  \right]
 + \rho C_{1\epsilon} f_1 \frac{\epsilon}{k} (P_k + C_{3\epsilon} G_b)
 - \rho C_{2\epsilon} f_2 \frac{\epsilon^2}{k}
 + \rho S_\epsilon,
\end{equation}

where

- $\nu = \mu/\rho$ and $\nu_t = \mu_t/\rho$ are the molecular and turbulent kinematic
  viscosities;
- $\rho$ is the density, which may have a functional dependency with temperature and pressure;
- $P_k$ is shear production of turbulent kinetic energy;
- $G_b$ is production (or destruction) of turbulence by buoyancy;
- $S_k$ and $S_\epsilon$ collect additional user- or model-specific source terms, including
  compressibility ($\gamma_M$, $k$ equation only), Yap ($\rho\gamma_Y$, $\epsilon$ equation, StandardTwoLayer and StandardLowRe only), and
  low-Re corrections ($G'$, $\epsilon$ equation);
- $f_1$ and $f_2$ are damping functions (typically $f_1 = 1$ in high-Re models, with
  $f_2$ used in low-Re variants).

The eddy viscosity is given by

\begin{equation}
\mu_t = \rho \, C_\mu^{\mathrm{eff}} \frac{k^2}{\epsilon},
\end{equation}

where $C_\mu^{\mathrm{eff}}$ can be:

- the classical constant $C_\mu = 0.09$ (standard model),
- a realizable value $C_\mu^\text{real}$ depending on $S^2$ and $W^2$,
- a two-layer or low-Re modified value, which depends on wall distance and Reynolds numbers
  specific to the chosen near-wall treatment.

Shear production is computed as

\begin{equation}
P_k = G_k = \mu_t S^2 - \frac{2}{3} \mu_t (\nabla \cdot \mathbf{u})^2
           - \frac{2}{3} \rho k \, \nabla \cdot \mathbf{u},
\end{equation}

with $S^2 = 2 S_{ij} S_{ij}$. The final two terms are optional compressibility
corrections that vanish when the flow is strictly incompressible. Buoyancy production is
modeled as

\begin{equation}
G_b = \beta(T) \, \frac{\mu_t}{\Pr_t} \, \mathbf{g} \cdot \nabla T,
\end{equation}

where $\Pr_t$ is the turbulent Prandtl number, $\mathbf{g}$ is gravity, $T$ is
temperature, and $\beta(T)$ is the thermal expansion coefficient. This follows from
$G_b = -\frac{\mu_t}{\rho \Pr_t} \mathbf{g} \cdot \nabla \rho$ with
$\nabla \rho = -\rho \beta \nabla T$: in a stable stratification ($\nabla T$ opposite to
$\mathbf{g}$) $G_b < 0$ and buoyancy destroys turbulence.

### Model constants (standard high-Re model)

The standard high-Re $k$–$\epsilon$ model uses the classical constants
[!cite](launder1974numerical):

- $C_\mu = 0.09$,
- $C_{1\epsilon} = 1.44$,
- $C_{2\epsilon} = 1.92$,
- $\sigma_k = 1.0$,
- $\sigma_\epsilon = 1.3$.

These can be overridden through the parameters of the relevant kernels
([kEpsilonTKESourceSink.md], [kEpsilonTKEDSourceSink.md], and [kEpsilonViscosity.md]) if a different
calibration is required. The realizable variants use the calibration of
[!cite](shih1995realizable), $C_{2} = 1.9$ and $\sigma_\epsilon = 1.2$ (see [#realizable-variant]).
The default of `C2_eps` in [kEpsilonTKEDSourceSink.md] follows the selected variant. $\sigma_k$ and
$\sigma_\epsilon$ are the `scaling_coeff` of the turbulent diffusion kernels (default 1) and must be
set in the input; they are not checked.

### Variant-specific forms

The implementation supports several model variants. For clarity we summarize the main
differences in the $\epsilon$ equation and $\mu_t$ definition.

#### Standard (high-Re)

- Eddy viscosity: $\mu_t = \rho C_\mu k^2 / \epsilon$.
- $f_1 = 1$, $f_2 = 1$.
- No explicit low-Re damping; near-wall behavior is typically handled by wall functions.
- Optional terms:

1. compressibility correction $\gamma_M$ (dilatation dissipation) as an extra sink in the $k$ equation,
2. nonlinear production $G_\text{nl}$ added to $P_k$.

#### StandardLowRe

- Same basic structure as the standard model, but:
- $\mu_t$ is multiplied by a low-Re damping function $f_\mu$,
- $f_2$ is given by `f2_SKE_LRe(Re_t)`,
- an extra production
  $G' = D f_2 (G_k + 2\mu k/d^2)\exp(-E\,Re_d^2)$ is added to $P_\epsilon$ when
  `use_low_re_Gprime = true` [!cite](lien1996low).
- These modifications allow integration to the wall ($y^+ \approx 1$) without wall
  functions, but require much finer near-wall meshes [!cite](wilcox2006turbulence).

#### StandardTwoLayer

- Near the wall, $\mu_t$ and the dissipation length scale $\ell_\epsilon$ are given by
  one of the two-layer prescriptions:
- `twoLayerWolfstein` [!cite](wolfstein1969velocity),
- `twoLayerNorrisReynolds` [!cite](norris1975oneequation).
- Away from the wall, the model reverts to the standard high-Re form.
- A blending function $\lambda(Re_d)$ [!cite](jongen1998simulation) gives a smooth transition
  for both $\mu_t$ and the $\epsilon$ equation (see [#two-layer]).

#### Realizable id=realizable-variant

- Uses the realizable eddy viscosity of [!cite](shih1995realizable):
  \begin{equation}
  C_\mu^\text{real} = \frac{1}{A_0 + A_s \, k U^* / \epsilon}, \qquad
  U^* = \sqrt{S_{ij}S_{ij} + \Omega_{ij}\Omega_{ij}},
  \end{equation}
  with $A_0 = 4.04$, $A_s = \sqrt{6}\cos\phi$, $\phi = \tfrac{1}{3}\arccos(\sqrt{6}\,W)$,
  $W = S_{ij}S_{jk}S_{ki}/\tilde{S}^3$ and $\tilde{S} = \sqrt{S_{ij}S_{ij}}$. The system-rotation
  term of $\Omega_{ij}$ is not included (stationary frame). In an equilibrium logarithmic layer
  $C_\mu^\text{real} \approx 0.09$.
- The $k$ equation is the same as in the standard model.
- The $\epsilon$ equation is replaced by
  \begin{equation}
  \rho C_1 S \epsilon - \rho C_{2} \frac{\epsilon^2}{k + \sqrt{\nu\epsilon}}
  + C_{1\epsilon}\frac{\epsilon}{k} C_{3\epsilon} G_b, \qquad
  C_1 = \max\left(0.43, \frac{\eta}{\eta+5}\right), \qquad \eta = S\frac{k}{\epsilon},
  \end{equation}
  with $S = \sqrt{2 S_{ij} S_{ij}}$. The production term does not contain $G_k$, and the destruction
  term uses `f2_RKE` $= k/(k+\sqrt{\nu\epsilon})$, so it stays finite as $k \to 0$.
- Model constants: $C_{1\epsilon} = 1.44$ (buoyancy term only), $C_{2} = 1.9$, $\sigma_k = 1.0$,
  $\sigma_\epsilon = 1.2$.

!alert warning title=Set $\sigma_\epsilon = 1.2$ for the realizable variants
$\sigma_\epsilon$ is not a parameter of [kEpsilonTKEDSourceSink.md]. It is the `scaling_coeff` of
the $\epsilon$ turbulent diffusion kernel (the one with `diffusion_coeff = mu_t`), or it is built
into a combined coefficient such as $\mu + \mu_t/\sigma_\epsilon$. Its default is 1 and nothing
checks it. `Realizable` and `RealizableTwoLayer` are calibrated with $\sigma_\epsilon = 1.2$ and
$C_2 = 1.9$ [!cite](shih1995realizable); the Standard variants use $\sigma_\epsilon = 1.3$.
`C2_eps` follows the variant automatically, but $\sigma_\epsilon$ has to be set in the input.

- Improves predictions for rotating, swirling, and separating flows while maintaining
  realizability of the normal stresses [!cite](shih1995realizable).

#### RealizableTwoLayer

- Combines the realizable $C_\mu$ formulation with a two-layer near-wall treatment.
- Often a good default for complex internal flows when mesh resolution allows a
  near-wall cluster with $1 \lesssim y^+ \lesssim 50$.

## Corrections implemented in the $k$–$\epsilon$ models

All correction functions are defined in [turbulenceModelingReference.md].
This section provides a summary of the corrections implemented in the $k$–$\epsilon$ models.

This section explains the available corrections, their motivation, and how they are used.

### Low-Reynolds-number corrections

- *Motivation*: Standard high-Re models with wall functions assume an equilibrium
  logarithmic layer and cannot accurately represent the viscous sublayer and buffer region.
  Low-Re formulations include damping functions that allow integration to the wall at
  $y^+ \approx 1$ [!cite](wilcox2006turbulence).
- *Implementation*:
  (a) $\mu_t \gets f_\mu \, \rho C_\mu k^2/\epsilon$ with `fmu_SKE_LRe`,
  (b) The $\epsilon$ equation uses $f_2 = f2_SKE_LRe$ or `f2_RKE`,
  (c) An additional production term $G'$ (via `computeGprime`) can be added to $S_\epsilon$.
- *Activation*:
  (a) `k_epsilon_variant = 'StandardLowRe'`,
  (b) `use_low_re_Gprime = true/false`,
  (c) near-wall mesh with $y^+ \approx 1$.

### Two-layer near-wall models id=two-layer

- *Motivation*: Two-layer models avoid the need to fully resolve the viscous sublayer
  while improving accuracy over simple wall functions. Near the wall, the $k$ equation is
  kept but $\epsilon$ and $\mu_t$ are given by algebraic length scales (a one-equation model);
  away from the wall, the standard or realizable $k$–$\epsilon$ equations are used
  [!cite](wolfstein1969velocity,norris1975oneequation). The two regions are blended
  smoothly following Jongen ([!cite](jongen1998simulation)).

#### Near-wall length scales

All near-wall quantities are functions of the wall-distance Reynolds number

\begin{equation}
Re_d = \frac{\sqrt{k}\, d}{\nu},
\end{equation}

where $d$ is the distance to the closest wall. The inner layer is described by a dissipation
length $\ell_\epsilon$ and a viscosity ratio $(\mu_t/\mu)_{2L}$, which give

\begin{equation}
\epsilon_{2L} = \frac{k^{3/2}}{\ell_\epsilon}, \qquad
\mu_{t,2L} = \left(\frac{\mu_t}{\mu}\right)_{2L} \mu.
\end{equation}

The available flavors are:

- Wolfstein ([!cite](wolfstein1969velocity)), with $c_\ell = 0.42\, C_\mu^{-3/4}$:
  \begin{equation}
  \ell_\epsilon = c_\ell\, d \left(1 - e^{-Re_d/(2 c_\ell)}\right), \qquad
  \left(\frac{\mu_t}{\mu}\right)_{2L} = 0.42\, C_\mu^{1/4}\, Re_d \left(1 - e^{-Re_d/70}\right),
  \end{equation}
- Norris–Reynolds ([!cite](norris1975oneequation)):
  \begin{equation}
  \ell_\epsilon = c_\ell\, d \frac{Re_d}{Re_d + 5.3}, \qquad
  \left(\frac{\mu_t}{\mu}\right)_{2L} = 0.42\, C_\mu^{1/4}\, Re_d \left(1 - e^{-Re_d/50.5}\right).
  \end{equation}

For both flavors the viscosity ratio is the one-equation
viscosity $\rho C_\mu \ell_\mu \sqrt{k} / \mu$, with $\ell_\mu = c_\ell\, d\, (1 - e^{-Re_d/A_\mu})$.

#### Blending function

The inner and outer layers are blended with the function of Jongen ([!cite](jongen1998simulation)):

\begin{equation}
\lambda = \tfrac12\left[1 + \tanh\!\left(\frac{Re_d - Re^*}{A}\right)\right],
\qquad A = \frac{\Delta Re}{\operatorname{atanh}(0.98)},
\end{equation}

so that $\lambda \to 0$ in the viscosity-affected region, $\lambda \to 1$ in the fully
turbulent region, and $\lambda$ goes from 0.01 to 0.99 over $Re^* \pm \Delta Re$.
$Re^*$ and $\Delta Re$ are the `Rey_star` (default 60) and `DeltaRey` (default 10) parameters of
[kEpsilonTwoLayerFunctorMaterial.md]. The same $\lambda$ is used for both the turbulent viscosity
and the $\epsilon$ equation, so the two never disagree about where the inner layer is.

#### Turbulent viscosity

\begin{equation}
\mu_t = \lambda\, \mu_{t,k\epsilon} + (1 - \lambda)\, \mu_{t,2L},
\end{equation}

where $\mu_{t,k\epsilon}$ is the standard or realizable outer-layer viscosity.

#### Dissipation equation

The $\epsilon$ equation is blended *at the discrete level*: each cell's row is a
$\lambda$-weighted combination of the transport equation and an algebraic relation that
drives $\epsilon$ toward $\epsilon_{2L}$:

\begin{equation}
\lambda \left[ \frac{\partial \rho \epsilon}{\partial t} + \nabla \cdot (\rho \mathbf{u} \epsilon)
- \nabla \cdot \left( \left(\mu + \frac{\mu_t}{\sigma_\epsilon}\right) \nabla \epsilon \right)
- P_\epsilon + D_\epsilon \right]
+ (1 - \lambda)\, a_{2L} \left( \epsilon - \epsilon_{2L} \right) = 0,
\qquad a_{2L} = \frac{\rho\, \epsilon_{2L}}{k}.
\end{equation}

The coefficient $a_{2L} = \rho \sqrt{k}/\ell_\epsilon$ is the inverse of the algebraic relaxation
time scale. It keeps the algebraic part of the row of the same magnitude as the destruction term
of the transport part. The limits are:

- $\lambda \to 1$: the standard $\epsilon$ transport equation,
- $\lambda \to 0$: $\epsilon = \epsilon_{2L}$ exactly, i.e. the classical two-layer model, where
  the $\epsilon$ equation is not solved in the inner layer.

The source and sink terms ($P_\epsilon$, $D_\epsilon$ and the algebraic part) are assembled by
[kEpsilonTKEDSourceSink.md]. The time-derivative, advection and diffusion terms are assembled by
the blended transport kernels [LinearFVTurbulentBlendedTimeDerivative.md],
[LinearFVTurbulentBlendedAdvection.md] and [LinearFVTurbulentBlendedDiffusion.md]. These scale
each row by the $\lambda$ of its cell.

In wall-bounded cells the transport rows are removed, as for the high-Re wall functions, and
$\epsilon$ is imposed as

\begin{equation}
\epsilon = \lambda\, \epsilon_\text{wall} + (1 - \lambda)\, \epsilon_{2L},
\end{equation}

where $\epsilon_\text{wall}$ is the wall-function value. Coarse meshes (first cell in the log
layer, $\lambda \to 1$) therefore recover the standard wall-function behavior, and meshes with
$y^+ \sim 1$ ($\lambda \to 0$) use the algebraic two-layer value.

The $k$ equation is *not* blended. It is solved in every cell with the regular turbulent
transport kernels and its usual destruction term $\rho \epsilon$. It sees the inner layer only
through the blended $\epsilon$ and $\mu_t$.

#### Activation

- `k_epsilon_variant = 'StandardTwoLayer'` or `'RealizableTwoLayer'` in [kEpsilonViscosity.md],
  [kEpsilonTKESourceSink.md] and [kEpsilonTKEDSourceSink.md],
- a `wall_distance` field (e.g. from [WallDistanceAux.md]),
- one [kEpsilonTwoLayerFunctorMaterial.md], which also selects
  `two_layer_flavor = 'Wolfstein' | 'NorrisReynolds'` and sets `Rey_star` and `DeltaRey`,
- `two_layer_lambda` and `two_layer_mu_t` set on [kEpsilonViscosity.md], and `two_layer_lambda`
  and `two_layer_epsilon` set on [kEpsilonTKEDSourceSink.md], all pointing to the outputs of that
  material,
- every transport kernel of $\epsilon$ replaced by its blended counterpart, with
  `blending_factor` set to the same `two_layer_lambda`, and the same `walls` list as
  [kEpsilonTKEDSourceSink.md] (this is not checked),
- the transport kernels of $k$ left as the regular, non-blended kernels with an empty `walls` list.

All of these rules except the `walls` match are checked when the simulation is set up, and an
inconsistent input stops with an error. A minimal setup is:

!listing test/tests/kEpsilon/special-cases/channel_ERCOFTAC_buoyant.i block=FunctorMaterials/two_layer

!listing test/tests/kEpsilon/special-cases/channel_ERCOFTAC_buoyant.i block=LinearFVKernels/TKED_advection

### Realizable formulation

- *Motivation*: The classical $k$–$\epsilon$ model can violate realizability constraints
  (e.g., negative normal stresses). The realizable formulation introduces a new expression
  for $C_\mu$ and a modified $\epsilon$ equation to enforce these constraints, improving
  predictions in rotating and separating flows [!cite](shih1995realizable).
- *Implementation*:
  (a) `Cmu_realizable` computes $C_\mu^\text{real}$ from the mean strain-rate and rotation-rate
    tensors, $k$ and $\epsilon$,
  (b) `C1_realizable` gives $C_1(\eta)$ for the production term $\rho C_1 S \epsilon$,
  (c) `f2_RKE` is always applied in the destruction term.
- *Activation*: `k_epsilon_variant = 'Realizable'` or `'RealizableTwoLayer'`.


### Nonlinear eddy-viscosity model

- *Motivation*: Linear eddy-viscosity models cannot reproduce some important features of
  anisotropic turbulence (secondary flows, normal stress differences) [!cite](pope2000turbulent).
  Nonlinear constitutive relations add higher-order dependencies on $\mathbf{S}$ and
  $\mathbf{W}$ [!cite](lien1996low).
- *Implementation*:
  (a) `computeTRANS_NL` constructs nonlinear tensor bases,
  (b) `computeGnl` computes additional production $G_\text{nl}$,
  (c) `Cmu_nonlinear` keeps $\mu_t$ consistent with the nonlinear stress model.
- *Activation*: `nonlinear_model = 'quadratic'` or `'cubic'` (in both the $k$ and $\epsilon$ kernels).
  Only available for the `Standard`, `StandardTwoLayer` and `StandardLowRe` variants.


### Yap correction

- *Motivation*: In separated and impinging flows, classical $k$–$\epsilon$ models may
  overpredict the turbulence length scale, leading to excessive eddy viscosity and overly
  diffused separation zones. Yap’s correction adds an extra source of $\epsilon$ to
  limit the growth of the length scale [!cite](yap1987correction).
- *Implementation*:
  (a) `computeGammaY` computes $\gamma_Y$ from $\ell = k^{3/2}/\epsilon$ and $\ell_\epsilon = C_l d$,
  (b) $\rho\gamma_Y$ is added to the $\epsilon$ equation in `kEpsilonTKEDSourceSink`.
- *Activation*: `use_yap = true` with `StandardTwoLayer` or `StandardLowRe` (other variants are
  an error) and a `wall_distance` functor.


### Buoyancy production

- *Motivation*: In natural or mixed convection, buoyancy can generate or suppress
  turbulence, especially in strongly stratified regions [!cite](wilcox2006turbulence).
- *Implementation*:
  (a) `computeGb` provides $G_b$ using gravity, temperature gradient, and
    expansion coefficient $\beta$,
  (b) $G_b$ is added to the $k$ equation and appears in the $\epsilon$ equation multiplied
    by $C_{3\epsilon}$.
- *Activation*:
  (a) `use_buoyancy = true`,
  (b) `gravity`, `temperature`, `beta`, and `Pr_t` set appropriately (in both the
    $k$ and $\epsilon$ kernels).


### Compressibility correction

- *Motivation*: In compressible turbulence, part of the dissipation is dilatational.
  Sarkar et al. model it as an extra dissipation proportional to the turbulent Mach number
  $M_t^2 = k/c^2$ [!cite](sarkar1991analysis).
- *Implementation*:
  (a) `computeGammaM` returns $\gamma_M = \rho C_M k \epsilon / c^2$, which is subtracted
    from the source of the $k$ equation only ($C_M = 2$, the default, recovers Sarkar's
    $2\rho\epsilon M_t^2$ with $M_t^2 = k/c^2$). It is applied after the $C_{PL}$ production
    limiter, so the limiter never removes it,
  (b) the same switch adds the dilatation terms
    $-\tfrac{2}{3}\rho k\nabla\cdot\mathbf{u} - \tfrac{2}{3}\mu_t(\nabla\cdot\mathbf{u})^2$ to $G_k$,
    which reach the $\epsilon$ equation through $C_{1\epsilon} G_k \epsilon/k$ (not in the
    realizable variants, whose $\epsilon$ production is $\rho C_1 S \epsilon$).
- *Activation*:
  (a) `use_compressibility = true` in **both** the $k$ and $\epsilon$ kernels,
  (b) a `speed_of_sound` functor in the $k$ kernel, consistent with the equation of state.


### Curvature / rotation correction

- *Motivation*: Simple eddy-viscosity models tend to mispredict turbulence levels in
  rotating channels, swirling jets, and flows over curved surfaces
  [!cite](spalart1997rotation). Arolla and Durbin extended the idea to scalar two-equation
  models with a factor $f_c$ that multiplies the production [!cite](arolla2013modeling).
- *Implementation*: `computeCurvatureFactor` evaluates
  \begin{equation}
  f_c = \min\left(C_\text{max},\;
  \frac{1}{C_{r1}\left(|\eta| - \eta\right) + \sqrt{1 - \min(C_{r2}\,\eta,\, 0.99)}}\right),
  \qquad
  \eta = \left(S_{ij}S_{ij} - W_{ij}W_{ij}\right) T^2,
  \end{equation}
  \begin{equation}
  T = \max(T_1, T_3), \qquad T_1 = \frac{k}{\epsilon}, \qquad
  T_2 = 6\sqrt{\frac{\nu}{\epsilon}}, \qquad T_3 = \left(T_1^{1.625}\, T_2\right)^{1/2.625},
  \end{equation}
  with $C_{r1} = 0.04645$, $C_{r2} = 0.25$ and $C_\text{max} = 1.25$.

  For the Realizable variants, $G_k$ in the $k$ equation and $\rho C_1 S\epsilon$ in the
  $\epsilon$ equation are multiplied by $f_c$. $f_c = 1$ in simple shear. Only the local
  rotation tensor is used; the streamline-curvature term of [!cite](wallin2002modelling)
  is not included.
- *Activation*: `use_curvature_correction = true` in both the $k$ and $\epsilon$ kernels
  (Realizable and RealizableTwoLayer only).


## Practical guidance for using each model

### Summary table

The table below provides a practical summary for choosing a $k$–$\epsilon$ variant in
OpenPronghorn.

| Variant              | Near-wall treatment                      | Typical first-cell $y^+$ | Mesh requirements                                    | Typical domain of application                                               | Example applications (OpenPronghorn context)                                   | Comments |
|----------------------|------------------------------------------|--------------------------|------------------------------------------------------|------------------------------------------------------------------------------|--------------------------------------------------------------------------------|----------|
| `Standard`           | High-Re, wall functions                  | $30 \lesssim y^+ \lesssim 300$ | Coarse first cell in log layer; few cells across BL  | Fully turbulent attached flows, mild pressure gradients                      | Straight channels, simple pipe flows, external bluff body with mild separation | Fast and robust; good initial choice for simple internal flows. |
| `StandardLowRe`      | Low-Re damping ($f_\mu$, $f_2$)         | $y^+ \approx 1$         | Strong wall clustering; many cells in viscous sublayer and buffer region | Heat-transfer-dominated internal flows, strong adverse pressure gradients | High-fidelity channel-flow benchmarks, heated rod bundles (resolved near wall) | More expensive; sensitive to near-wall mesh quality and wall-distance accuracy. |
| `StandardTwoLayer`   | Two-layer (Wolfstein/NorrisReynolds)     | $1 \lesssim y^+ \lesssim 100$ | Moderate wall clustering; smooth near-wall blending  | Industrial internal flows, moderate separation, HVAC and ductwork           | Reactor auxiliary systems, plenum mixing with moderate separation             | Good compromise between cost and accuracy; pick `two_layer_flavor` based on target data. |
| `Realizable`         | High-Re, wall functions                  | $30 \lesssim y^+ \lesssim 300$ | Similar to `Standard`                                | Rotating, swirling, and separating flows                                     | Mixing tees, swirling jets, complex manifolds with recirculation              | Often preferred to `Standard` when rotation or separation is important. |
| `RealizableTwoLayer` | Realizable + two-layer near-wall model   | $1 \lesssim y^+ \lesssim 100$ | Moderate wall clustering plus two-layer blending     | General-purpose model for complex internal flows with heat transfer          | Reactor primary loops with moderate mesh resolution near walls                | Recommended default in many industrial-like configurations when mesh allows $y^+ \lesssim 50$. |

### Model setup id=model-setup

A $k$–$\epsilon$ model is set up across several objects: [kEpsilonTKESourceSink.md] ($k$),
[kEpsilonTKEDSourceSink.md] ($\epsilon$), [kEpsilonViscosity.md] ($\mu_t$), and, for the two-layer
variants, [kEpsilonTwoLayerFunctorMaterial.md]. The table shows which options each variant accepts.

| Option | `Standard` | `StandardLowRe` | `StandardTwoLayer` | `Realizable` | `RealizableTwoLayer` |
|---|---|---|---|---|---|
| `nonlinear_model` ($k$, $\epsilon$) | yes | yes | yes | error | error |
| `use_yap` ($\epsilon$) | error | yes | yes | error | error |
| `use_low_re_Gprime` ($\epsilon$) | error | yes | error | error | error |
| `curvature_model` ($k$, $\epsilon$) | error | error | error | yes | yes |
| `use_buoyancy` ($k$, $\epsilon$) | yes | yes | yes | yes | yes |
| `use_compressibility` ($k$, $\epsilon$) | yes | yes | yes | yes (1) | yes (1) |
| `Cd0`, `Cd1`, `Cd2` ($\mu_t$) | error | yes | error | error | error |
| `A0` ($\mu_t$) | error | error | error | yes | yes |
| `two_layer_lambda`, `two_layer_epsilon`, `two_layer_mu_t` | error | error | required | error | required |
| `wall_distance` ($\epsilon$, $\mu_t$) | -- | required | required ($\epsilon$) | -- | required ($\epsilon$) |
| `walls` ($k$) | optional | optional | required | optional | required |

"error" means the input stops with an error message. (1) In the realizable variants the
$\epsilon$ production contains no $G_k$, so the flag only acts in the $k$ equation.

Some constants are accepted but have no effect outside the model they belong to:
`C_lowRe` (StandardLowRe only), `D_lowRe` and `E_lowRe` (only with `use_low_re_Gprime`), `Cw`
(only with `use_yap`), `C3_eps`, `Pr_t` and `gravity` in [kEpsilonTKEDSourceSink.md] (only with
`use_buoyancy`), and `Ct` in [kEpsilonViscosity.md] (only with `scale_limiter = standard`).

!alert warning title=Keep the shared settings identical
The following settings are read separately by each object, and nothing checks that they agree:
`k_epsilon_variant`, `nonlinear_model`, `use_compressibility`, `use_buoyancy`, `curvature_model`,
`gravity`, `Pr_t`, `temperature`, `beta`, `C_mu`, `walls` and `wall_treatment`. They must be the
same in [kEpsilonTKESourceSink.md], [kEpsilonTKEDSourceSink.md] and [kEpsilonViscosity.md], where
they apply. `C_mu` must also match in [kEpsilonTwoLayerFunctorMaterial.md] and in the wall-function
boundary condition of $\mu_t$, and `walls` in the $\epsilon$ transport kernels. A mismatch, for
example a different `k_epsilon_variant` in the $k$ and $\epsilon$ kernels, runs without error but
does not solve any of the models described here. The defaults of `wall_treatment` also differ
(`neq` in the $k$ and $\epsilon$ kernels and the wall-function boundary condition, `eq_newton` in
[kEpsilonViscosity.md]), so set it explicitly in every object. Defining each shared value once
as an input-file variable, as the test inputs do, avoids mismatches.

### Additional configuration guidelines

1. *Wall-distance field*

- Use `WallDistanceAux` to compute the distance to solid walls.
- Ensure that the `walls` list is consistent across `WallDistanceAux`,
     `kEpsilonTKESourceSink`, `kEpsilonTKEDSourceSink`, and `kEpsilonViscosity`.

2. *Coupling $k$, $\epsilon$, and $\mu_t$*

- Always include all three components: `kEpsilonTKESourceSink` (for $k$),
     `kEpsilonTKEDSourceSink` (for $\epsilon$), and `kEpsilonViscosity` (for $\mu_t$).
- The `k_epsilon_variant` parameter, and the other shared settings listed in
     [#model-setup], must be the same in all of them.

3. *Incompressible vs compressible*

- For incompressible or low-Mach flows, set `use_compressibility = false` and omit
     `speed_of_sound`.
- For compressible flows, provide a consistent `speed_of_sound` functor (to the $k$ kernel) and consider
     tighter nonlinear and linear solver tolerances.

4. *Buoyancy-driven and stratified flows*

- Pay attention to `Pr_t`, `beta`, and the reference temperature/state.
- Post-process $P_k$, $G_b$, and $S_\epsilon$ terms where possible to verify the
     qualitative behavior of the model.

5. *Enabling corrections judiciously*

- Nonlinear and curvature corrections add complexity and stiffness. Use them when
     there is clear physical motivation (strong curvature, swirl, or secondary flows) and
     when validation data are available.
- Yap and low-Re corrections are valuable in separated flows and near-wall heat
     transfer, but may not be necessary in simple attached turbulence.

6. *Validation and calibration*

- The default constants are consistent with widely used industrial models
     [!cite](launder1974numerical,shih1995realizable).
- For critical design calculations, validate against experiments or high-fidelity
     simulations (DNS/LES) [!cite](pope2000turbulent,wilcox2006turbulence).
- If any model coefficients are retuned, document these changes and the validation
     cases used.
