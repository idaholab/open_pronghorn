# $k$–$\epsilon$ Turbulence Modeling in OpenPronghorn

This page documents the Reynolds-averaged Navier–Stokes (RANS) turbulence-modeling
capabilities implemented in OpenPronghorn. It focuses on the $k$–$\epsilon$ family of
eddy-viscosity models and practical guidance for choosing and using each model in applications.

For general background on turbulence and RANS closures, see Pope [!cite](pope2000turbulent)
and Wilcox [!cite](wilcox2006turbulence). The baseline $k$–$\epsilon$ model follows
Launder and Spalding [!cite](launder1974numerical), with additional options inspired by
realizable, low-Reynolds-number, two-layer, curvature/rotation, buoyancy, and
compressibility corrections [!cite](shih1995realizable,spalart1997rotation,yap1987correction,wolfstein1969velocity,norris1975oneequation).


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
 + \rho P_k + \rho G_b - \rho \epsilon + \rho S_k,
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
  compressibility ($\gamma_M$), Yap ($\gamma_Y$), and low-Re corrections ($G'$);
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
G_b = - \frac{\mu_t}{\Pr_t} \, \mathbf{g} \cdot \nabla T \, \beta(T),
\end{equation}

where $\Pr_t$ is the turbulent Prandtl number, $\mathbf{g}$ is gravity, $T$ is
temperature, and $\beta(T)$ is the thermal expansion coefficient.

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
calibration is required.

### Variant-specific forms

The implementation supports several model variants. For clarity we summarize the main
differences in the $\epsilon$ equation and $\mu_t$ definition.

#### Standard (high-Re)

- Eddy viscosity: $\mu_t = \rho C_\mu k^2 / \epsilon$.
- $f_1 = 1$, $f_2 = 1$.
- No explicit low-Re damping; near-wall behavior is typically handled by wall functions.
- Optional terms:

1. compressibility correction $\gamma_M$ in $S_\epsilon$,
2. Yap correction $\gamma_Y$ in $S_\epsilon$,
3. curvature factor $f_c$ multiplying $P_k$,
4. nonlinear production $G_\text{nl}$ added to $P_k$.

#### StandardLowRe

- Same basic structure as the standard model, but:
- $\mu_t$ is multiplied by a low-Re damping function $f_\mu$,
- $f_2$ is given by `f2_SKE_LRe(Re_t)`,
- an additional $G'$ term may appear in the $\epsilon$ equation, scaled by
    a low-Re damping factor.
- These modifications allow integration to the wall ($y^+ \approx 1$) without wall
  functions, but require much finer near-wall meshes [!cite](wilcox2006turbulence).

#### StandardTwoLayer

- Near the wall, $\mu_t$ and the dissipation length scale $\ell_\epsilon$ are given by
  one of the two-layer prescriptions:
- `twoLayerWolfstein` [!cite](wolfstein1969velocity),
- `twoLayerNorrisReynolds` [!cite](norris1975oneequation),
- `twoLayerXu` [!cite](xu1998new).
- Away from the wall, the model reverts to the standard high-Re form.
- A blending function based on wall-distance Reynolds numbers provides a smooth transition.

#### Realizable

- Uses the realizable form of $C_\mu$ [!cite](shih1995realizable):
  \begin{equation}
  C_\mu^\text{real} =
  \frac{C_{a0}}{C_{a1} + C_{a2} \bar{S} + C_{a3} \bar{W}},
  \end{equation}
  where $\bar{S} = (k/\epsilon) \sqrt{S^2}$ and $\bar{W} = (k/\epsilon) \sqrt{W^2}$.
- Optionally uses a realizable low-Re damping function `f2_RKE`.
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
  (a) $\mu_t \gets f_\mu \, \rho k^2/\epsilon$ with `fmu_SKE_LRe` (standard) or
    combined realizable damping,
  (b) The $\epsilon$ equation uses $f_2 = f2_SKE_LRe$ or `f2_RKE`,
  (c) An additional production term $G'$ (via `computeGprime`) can be added to $S_\epsilon$.
- *Activation*:
  (a) `k_epsilon_variant = 'StandardLowRe'`,
  (b) `use_low_re_Gprime = true/false`,
  (c) near-wall mesh with $y^+ \approx 1$.

### Two-layer near-wall models

- *Motivation*: Two-layer models avoid the need to fully resolve the viscous sublayer
  while improving accuracy over simple wall functions. Near the wall, $\epsilon$ and
  $\mu_t$ are driven by empirical formulas; away from the wall, the standard or realizable
  $k$–$\epsilon$ equations are used
  [!cite](wolfstein1969velocity,norris1975oneequation).
- *Implementation*:
  (a) The wall distance $d$ is obtained from `WallDistanceAux`,
  (b) `twoLayerWolfstein`, `twoLayerNorrisReynolds`, or `twoLayerXu` compute $\ell_\epsilon$
    and $\mu_t/\mu$ in the near-wall region,
  (c) A blending function transitions between near-wall and outer-layer behavior.
- *Activation*:
  (a) `k_epsilon_variant = 'StandardTwoLayer'` or `'RealizableTwoLayer'`,
  (b) `two_layer_flavor = 'Wolfstein' | 'NorrisReynolds' | 'Xu'`,
  (c) `wall_distance` auxiliary variable and consistent `walls` list.

### Realizable formulation

- *Motivation*: The classical $k$–$\epsilon$ model can violate realizability constraints
  (e.g., negative normal stresses). The realizable formulation introduces a new expression
  for $C_\mu$ and a modified $\epsilon$ equation to enforce these constraints, improving
  predictions in rotating and separating flows [!cite](shih1995realizable).
- *Implementation*:
  (a) `Cmu_realizable` provides $C_\mu^\text{real}$ based on $S^2$, $W^2$, $k$, and
    $\epsilon$.
  (b) Optional low-Re damping via `f2_RKE` may be applied.
- *Activation*: `k_epsilon_variant = 'Realizable'` or `'RealizableTwoLayer'`.


### Nonlinear eddy-viscosity model

- *Motivation*: Linear eddy-viscosity models cannot reproduce some important features of
  anisotropic turbulence (secondary flows, normal stress differences) [!cite](pope2000turbulent).
  Nonlinear constitutive relations add higher-order dependencies on $\mathbf{S}$ and
  $\mathbf{W}$.
- *Implementation*:
  (a) `computeTRANS_NL` constructs nonlinear tensor bases,
  (b) `computeGnl` computes additional production $G_\text{nl}$,
  (c) `Cmu_nonlinear` keeps $\mu_t$ consistent with the nonlinear stress model.
- *Activation*: `nonlinear_model = 'quadratic'` or `'cubic'` (in the $k$–$\epsilon$ kernels).


### Yap correction

- *Motivation*: In separated and impinging flows, classical $k$–$\epsilon$ models may
  overpredict the turbulence length scale, leading to excessive eddy viscosity and overly
  diffused separation zones. Yap’s correction adds an extra sink term for $\epsilon$ to
  limit the growth of the length scale [!cite](yap1987correction).
- *Implementation*:
  (a) `computeGammaY` computes $\gamma_Y$ based on $\ell$, $\ell_\epsilon$, and model
    constants,
  (b) $\gamma_Y$ is added to $S_\epsilon$ in `kEpsilonTKEDSourceSink` when enabled.
- *Activation*:
  (a) `use_yap = true`,
  (b) consistent wall-distance and two-layer settings if applicable.


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
  (b) `gravity`, `T_fluid`, `beta`, and `Pr_t` set appropriately.


### Compressibility correction

- *Motivation*: At high Mach numbers, dilatational effects can alter the turbulence
  structure and time scales. Compressibility corrections add an extra term to the
  $\epsilon$ equation to account for these effects [!cite](wilcox2006turbulence).
- *Implementation*:
  (a) `computeGammaM` computes $\gamma_M$ based on $k$, $\epsilon$, and local speed of
    sound,
  (b) $\gamma_M$ is added to $S_\epsilon$ when `use_compressibility = true`.
- *Activation*:
  (a) `use_compressibility = true`,
  (b) a `speed_of_sound` functor is provided and consistent with the equation of state.


### Curvature / rotation correction

- *Motivation*: Simple eddy-viscosity models tend to mispredict turbulence levels in
  rotating channels, swirling jets, and flows over curved surfaces. Spalart and Shur
  proposed a correction factor $f_c$ that modifies production based on curvature/rotation
  [!cite](spalart1997rotation).
- *Implementation*:
  (a) `computeCurvatureFactor` uses $S^2$ and $W^2$ to build $f_c$,
  (b) $P_k$ is replaced by $f_c P_k$ in both the $k$ and $\epsilon$ equations when
    curvature corrections are enabled.
- *Activation*:
  (a) `use_curvature_correction = true`,
  (b) `curvature_model = 'standard'` (or similar option as implemented).


## Practical guidance for using each model

### Summary table

The table below provides a practical summary for choosing a $k$–$\epsilon$ variant in
OpenPronghorn.

| Variant              | Near-wall treatment                      | Typical first-cell $y^+$ | Mesh requirements                                    | Typical domain of application                                               | Example applications (OpenPronghorn context)                                   | Comments |
|----------------------|------------------------------------------|--------------------------|------------------------------------------------------|------------------------------------------------------------------------------|--------------------------------------------------------------------------------|----------|
| `Standard`           | High-Re, wall functions                  | $30 \lesssim y^+ \lesssim 300$ | Coarse first cell in log layer; few cells across BL  | Fully turbulent attached flows, mild pressure gradients                      | Straight channels, simple pipe flows, external bluff body with mild separation | Fast and robust; good initial choice for simple internal flows. |
| `StandardLowRe`      | Low-Re damping ($f_\mu$, $f_2$)         | $y^+ \approx 1$         | Strong wall clustering; many cells in viscous sublayer and buffer region | Heat-transfer-dominated internal flows, strong adverse pressure gradients | High-fidelity channel-flow benchmarks, heated rod bundles (resolved near wall) | More expensive; sensitive to near-wall mesh quality and wall-distance accuracy. |
| `StandardTwoLayer`   | Two-layer (Wolfstein/NorrisReynolds/Xu)  | $1 \lesssim y^+ \lesssim 100$ | Moderate wall clustering; smooth near-wall blending  | Industrial internal flows, moderate separation, HVAC and ductwork           | Reactor auxiliary systems, plenum mixing with moderate separation             | Good compromise between cost and accuracy; pick `two_layer_flavor` based on target data. |
| `Realizable`         | High-Re, wall functions                  | $30 \lesssim y^+ \lesssim 300$ | Similar to `Standard`                                | Rotating, swirling, and separating flows                                     | Mixing tees, swirling jets, complex manifolds with recirculation              | Often preferred to `Standard` when rotation or separation is important. |
| `RealizableTwoLayer` | Realizable + two-layer near-wall model   | $1 \lesssim y^+ \lesssim 100$ | Moderate wall clustering plus two-layer blending     | General-purpose model for complex internal flows with heat transfer          | Reactor primary loops with moderate mesh resolution near walls                | Recommended default in many industrial-like configurations when mesh allows $y^+ \lesssim 50$. |

### Additional configuration guidelines

1. *Wall-distance field*

- Use `WallDistanceAux` to compute the distance to solid walls.
- Ensure that the `walls` list is consistent across `WallDistanceAux`,
     `kEpsilonTKESourceSink`, `kEpsilonTKEDSourceSink`, and `kEpsilonViscosity`.

2. *Coupling $k$, $\epsilon$, and $\mu_t$*

- Always include all three components: `kEpsilonTKESourceSink` (for $k$),
     `kEpsilonTKEDSourceSink` (for $\epsilon$), and `kEpsilonViscosity` (for $\mu_t$).
- The `k_epsilon_variant` parameter should be the same in all of them.

3. *Incompressible vs compressible*

- For incompressible or low-Mach flows, set `use_compressibility = false` and omit
     `speed_of_sound`.
- For compressible flows, provide a consistent `speed_of_sound` functor and consider
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
