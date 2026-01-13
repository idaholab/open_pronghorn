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


# TurbulenceMethods: complete list of utilities (equation-level detail)

This document focuses exclusively on the utilities implemented in
the code and explains, for each function, exactly which mathematical
expressions are implemented and how they are built from the input quantities.


## 1. `computeStrainRotationInvariants`

This function computes *velocity-gradient-based invariants* used in turbulence modeling:
the strain-rate invariant $S^2$, the rotation-rate invariant $W^2$, and the divergence
$\nabla \cdot \mathbf{u}$.

1. It first obtains the gradient of the streamwise component $u$:
   \begin{equation}
   \nabla u = \frac{\partial u}{\partial x_j}
            = \left( \frac{\partial u}{\partial x},
                      \frac{\partial u}{\partial y},
                      \frac{\partial u}{\partial z} \right).
   \end{equation}
   Analogously, when provided, it obtains $\nabla v$ and $\nabla w$ for the other
   velocity components.

2. The code then defines the *velocity-gradient components* as
   \begin{equation}
   \nabla u = \frac{\partial u_i}{\partial x_j}
   \begin{aligned}
   &\text{du}_x = \frac{\partial u}{\partial x}, \quad
    \text{du}_y = \frac{\partial u}{\partial y}, \quad
    \text{du}_z = \frac{\partial u}{\partial z}, \\
   &\text{dv}_x = \frac{\partial v}{\partial x}, \quad
    \text{dv}_y = \frac{\partial v}{\partial y}, \quad
    \text{dv}_z = \frac{\partial v}{\partial z}, \\
   &\text{dw}_x = \frac{\partial w}{\partial x}, \quad
    \text{dw}_y = \frac{\partial w}{\partial y}, \quad
    \text{dw}_z = \frac{\partial w}{\partial z}.
   \end{aligned}
   \end{equation}

3. The *divergence* of the velocity is computed as
   \begin{equation}
   \nabla \cdot \mathbf{u} = \frac{\partial u_i}{\partial x_j}
                           = \frac{\partial u}{\partial x}
                           + \frac{\partial v}{\partial y}
                           + \frac{\partial w}{\partial z},
   \end{equation}
   with terms omitted in lower dimensions. In code this is
   `div_u = dux + dvy + dwz` where the absent components default to zero.

4. The *symmetric strain-rate tensor* $S_{ij} = \frac{1}{2} \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right)$ is constructed as
   \begin{equation}
   S_{xx} = \frac{\partial u}{\partial x}, \quad
   S_{yy} = \frac{\partial v}{\partial y}, \quad
   S_{zz} = \frac{\partial w}{\partial z},
   \end{equation}
   \begin{equation}
   S_{xy} = S_{yx} = \frac{1}{2}
   \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right),
   \end{equation}
   \begin{equation}
   S_{xz} = S_{zx} = \frac{1}{2}
   \left( \frac{\partial u}{\partial z} + \frac{\partial w}{\partial x} \right),
   \end{equation}
   \begin{equation}
   S_{yz} = S_{zy} = \frac{1}{2}
   \left( \frac{\partial v}{\partial z} + \frac{\partial w}{\partial y} \right).
   \end{equation}

5. The *antisymmetric rotation tensor* $W_{ij} = \frac{1}{2} \left( \frac{\partial u_i}{\partial x_j} - \frac{\partial u_j}{\partial x_i} \right)$ is
   \begin{equation}
   W_{xy} = -W_{yx} = \frac{1}{2}
   \left( \frac{\partial u}{\partial y} - \frac{\partial v}{\partial x} \right),
   \end{equation}
   \begin{equation}
   W_{xz} = -W_{zx} = \frac{1}{2}
   \left( \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x} \right),
   \end{equation}
   \begin{equation}
   W_{yz} = -W_{zy} = \frac{1}{2}
   \left( \frac{\partial v}{\partial z} - \frac{\partial w}{\partial y} \right).
   \end{equation}

6. Finally the *invariants* returned are
   \begin{equation}
   S^2 = 2 S_{ij} S_{ij} =
   2\left(S_{xx}^2 + S_{yy}^2 + S_{zz}^2
         + 2(S_{xy}^2 + S_{xz}^2 + S_{yz}^2)\right),
   \end{equation}
   \begin{equation}
   W^2 = 2 W_{ij} W_{ij} =
   4(W_{xy}^2 + W_{xz}^2 + W_{yz}^2).
   \end{equation}

The struct `StrainRotationInvariants` stores `S2 = S^2`, `W2 = W^2`, and `div_u` for use by
other turbulence utilities.


## 2. `cl_from_Cmu`

This helper computes a *near-wall length-scale constant* $c_\ell$ as a function of
$C_\mu$:
\begin{equation}
c_\ell(C_\mu) = 0.42 \, C_\mu^{-0.75}.
\end{equation}

This constant is used by the two-layer models (`twoLayerWolfstein`, `twoLayerNorrisReynolds`).


## 3. `twoLayerWolfstein`

This function implements a *Wolfstein-type two-layer near-wall model*. It returns a
`TwoLayerLengths` struct with

- `l_eps`: dissipation length scale $\ell_\epsilon$,
- `mu_ratio`: viscosity ratio $\mu_t / \mu$.

Using $c_\ell = c_\ell(C_\mu)$ from `cl_from_Cmu`,
\begin{equation}
\ell_\epsilon = c_\ell \, d \left(1 - e^{- \mathrm{Re}_d / (2 c_\ell)}\right),
\end{equation}
\begin{equation}
\frac{\mu_t}{\mu} = 0.42 \, \mathrm{Re}_d^{0.25}
\left( 1 - e^{-\mathrm{Re}_d / 70} \right),
\end{equation}
where

- $d$ is the wall distance,
- $\mathrm{Re}_d$ is a Reynolds number based on wall distance.


## 4. `twoLayerNorrisReynolds`

This function implements a *Norris–Reynolds-type two-layer model*. It returns

\begin{equation}
\ell_\epsilon = c_\ell \, d \frac{\mathrm{Re}_d}{\mathrm{Re}_d + 5.3},
\end{equation}
\begin{equation}
\frac{\mu_t}{\mu} = 0.42 \, \mathrm{Re}_d^{0.25}
\left( 1 - e^{-\mathrm{Re}_d / 50.5} \right),
\end{equation}

with $c_\ell = c_\ell(C_\mu)$ as before.


## 5. `twoLayerXu`

This function implements an *Xu-type two-layer correlation* using a modified wall
coordinate $y_\nu^\ast$ (denoted `yv_star` in the code). It returns

\begin{equation}
\ell_\epsilon = \frac{8.8 \, d}
                      {1 + \frac{10}{y_\nu^\ast} + 5.15 \times 10^{-2} y_\nu^\ast},
\end{equation}
\begin{equation}
\frac{\mu_t}{\mu} = \frac{0.544 \, y_\nu^\ast}
                         {1 + 5.025 \times 10^{-4} \, y_\nu^{\ast\,1.65}}.
\end{equation}

Note that `Cmu` and `Re_d` are not used in this formulation.


## 6. `f2_SKE_LRe`

This function provides a *low-Reynolds-number damping function* $f_2$ for the
standard $k\text{–}\epsilon$ model:
\begin{equation}
f_2(\mathrm{Re}_t) = 1 - C \exp(-\mathrm{Re}_t^2),
\end{equation}
where $C$ is a model constant and $\mathrm{Re}_t$ is a turbulent Reynolds number
(typically $k^2/(\nu \epsilon)$).


## 7. `fmu_SKE_LRe`

This function computes a *low-Re damping factor* $f_\mu$ for the eddy viscosity:
\begin{equation}
\begin{aligned}
\text{arg} &= C_{d0} \sqrt{\mathrm{Re}_d}
           + C_{d1} \mathrm{Re}_d
           + C_{d2} \mathrm{Re}_d^2, \\\\
f_\mu(\mathrm{Re}_d) &= 1 - e^{-\text{arg}}.
\end{aligned}
\end{equation}

This is used to modify $\mu_t$ as $\mu_t \gets f_\mu \mu_t$ for the first cell near the wall.


## 8. `f2_RKE`

This is another *damping function* $f_2$ (used in realizable or low-Re variants):
\begin{equation}
f_2(k,\nu,\epsilon) = \frac{k}{k + \sqrt{\nu \epsilon}}.
\end{equation}

It approaches 0 near the wall (small $k$), and 1 away from the wall.


## 9. `Cmu_realizable`

This function computes a *realizable* eddy-viscosity constant $C_\mu$ based on the
invariants $S^2$ and $W^2$:

1. First compute the ratio
   \begin{equation}
   \frac{k}{\epsilon} = k_\epsilon,
   \end{equation}
   and then
   \begin{equation}
   \bar{S} = k_\epsilon \sqrt{S^2}, \qquad
   \bar{W} = k_\epsilon \sqrt{W^2}.
   \end{equation}

2. The realizable $C_\mu$ is then
   \begin{equation}
   C_\mu^\text{real} =
   \frac{C_{a0}}{C_{a1} + C_{a2} \bar{S} + C_{a3} \bar{W}}.
   \end{equation}

with the constants pre-defined realizability constants
$C_{a0} = 0.667$, $C_{a1} = 1.25$, $C_{a2} = 1.0$, and $C_{a3} = 0.9$.

This quantity is used in realizable $k\text{–}\epsilon$ models to define $\mu_t$.


## 10. `Cmu_nonlinear`

This function provides a *nonlinear-model-compatible* $C_\mu$. If $k \le 0$ or
$\epsilon \le 0$, it returns 0. Otherwise:

1. Use `inv.S2` and `inv.W2` from `computeStrainRotationInvariants` and clamp them to
   non-negative values.
2. Compute
   \begin{equation}
   \bar{S} = \frac{k}{\epsilon} \sqrt{\max(S^2, 0)}, \qquad
   \bar{W} = \frac{k}{\epsilon} \sqrt{\max(W^2, 0)},
   \end{equation}
3. Then
   \begin{equation}
   C_\mu^\text{NL} =
   \frac{C_{a0}}{C_{a1} + C_{a2} \bar{S} + C_{a3} \bar{W}}.
   \end{equation}

This has the same functional form as `Cmu_realizable`, but uses the invariants packaged in
`StrainRotationInvariants` and is used with nonlinear constitutive relations.


## 11. `computeGk`

This function implements the *shear-production term* $G_k$ for the $k$-equation.

- Always includes the basic term
  \begin{equation}
  G_k = \mu_t S^2.
  \end{equation}

- If `include_compressibility_terms == true`, it further subtracts the compressibility
  corrections:
  \begin{equation}
  G_k \;\gets\; G_k
    - \frac{2}{3}\left( \rho k \, \nabla \cdot \mathbf{u}
                     + \mu_t (\nabla \cdot \mathbf{u})^2 \right).
  \end{equation}

So the full expression is
\begin{equation}
G_k = \mu_t S^2
      - \frac{2}{3}
      \Bigl( \rho k \, \nabla \cdot \mathbf{u}
           + \mu_t (\nabla \cdot \mathbf{u})^2 \Bigr)
\quad\text{(if compressibility included)}.
\end{equation}

Note, however, that for incompressible flows $\nabla \cdot \mathbf{u} = 0$.
So, `include_compressibility_terms` should have no effect on the shear-production term.


## 12. `computeGb`

This function computes the *buoyancy production* term $G_b$ in the $k$-equation as
\begin{equation}
G_b = \beta \, \frac{\mu_t}{\Pr_t} \, (\nabla T \cdot \mathbf{g}),
\end{equation}
where

- $\beta$ is the thermal expansion coefficient,
- $\Pr_t$ is the turbulent Prandtl number,
- $\mathbf{g}$ is the gravity vector,
- $\nabla T$ is the temperature gradient.


## 13. `computeGammaM`

This function returns the *compressibility correction* $\gamma_M$ added to the
$\epsilon$-equation:

\begin{equation}
\gamma_M = \rho \, C_M \, \frac{k \, \epsilon}{c^2},
\end{equation}
where

- $C_M$ is a model constant,
- $c$ is the speed of sound.


## 14. `computeGammaY`

This function implements the *Yap correction* $\gamma_Y$ for the $\epsilon$
equation.

1. If $k \le 0$, return 0.
2. Otherwise compute
   \begin{equation}
   r = \frac{\ell}{\ell_\epsilon},
   \quad
   \text{tmp} = (r - 1) r^2,
   \end{equation}
   where $\ell = \frac{k^{3/2}}{\epsilon}$ is the bulk eddy lenght scale, then take
   \begin{equation}
   \text{val} = \max(\text{tmp}, 0).
   \end{equation}

3. The final correction is
   \begin{equation}
   \gamma_Y = C_w \, \frac{\epsilon^2}{k} \, \text{val}.
   \end{equation}

So $\gamma_Y$ injects extra dissipation whenever $\ell / \ell_\epsilon > 1$.


## 15. `computeGprime`

This function computes an additional *low-Re production* term $G'$ for the
$\epsilon$-equation.

1. First build
   \begin{equation}
   \text{term} = G_k + 2 \mu_t \frac{k}{d^2},
   \end{equation}
   where

   - $G_k$ is the shear-production term,
   - $d$ is wall distance.

2. Then
   \begin{equation}
   G' = D \, f_2 \, \text{term} \, \exp(-E \, \mathrm{Re}_d^2),
   \end{equation}
   where $D$ and $E$ are constants, $f_2$ is typically a low-Re damping function,
   and $\mathrm{Re}_d$ is a Reynolds number based on wall distance.


## 16. `computeVelocityGradient`

This function constructs the full *velocity gradient tensor*
$\nabla \mathbf{u} = \partial u_i/\partial x_j$ as a `RankTwoTensor`.

- It queries the gradients of $u$, and, when present, $v$ and $w$ and fills:
  \begin{equation}
  (\nabla \mathbf{u})_{00} = \frac{\partial u}{\partial x}, \quad
  (\nabla \mathbf{u})_{01} = \frac{\partial u}{\partial y}, \quad
  (\nabla \mathbf{u})_{02} = \frac{\partial u}{\partial z},
  \end{equation}
  \begin{equation}
  (\nabla \mathbf{u})_{10} = \frac{\partial v}{\partial x}, \quad
  (\nabla \mathbf{u})_{11} = \frac{\partial v}{\partial y}, \quad
  (\nabla \mathbf{u})_{12} = \frac{\partial v}{\partial z},
  \end{equation}
  \begin{equation}
  (\nabla \mathbf{u})_{20} = \frac{\partial w}{\partial x}, \quad
  (\nabla \mathbf{u})_{21} = \frac{\partial w}{\partial y}, \quad
  (\nabla \mathbf{u})_{22} = \frac{\partial w}{\partial z}.
  \end{equation}

This tensor is then used by the nonlinear constitutive relation and curvature correction.


## 17. Nonlinear constitutive relation utilities

The remaining utilities implement *quadratic and cubic nonlinear eddy-viscosity models*
and associated terms.

There are fixed model coefficients in an anonymous namespace:

- Quadratic/cubic model coefficients: `CNL1`–`CNL7`,
- `Ca0`–`Ca3` for the nonlinear $C_\mu$ expression. These are used internally by
  `computeTRANS_NL` and `computeGnl`.


### 17.1 `computeTRANS_NL`

This function builds the *nonlinear part of the Reynolds-stress tensor* associated with
a quadratic or cubic constitutive relation.

1. If `model == None` or any of $\mu_t, k, \epsilon$ are non-positive, it returns
   the zero tensor.

2. Construct the *strain* and *rotation* tensors from the velocity gradient:
   \begin{equation}
   S = \frac{1}{2} \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right), \qquad
   W = \frac{1}{2} \left( \nabla \mathbf{u} - (\nabla \mathbf{u})^T \right).
   \end{equation}

3. Use the invariants in `inv` to compute
   \begin{equation}
   S\cdot S = \frac{1}{2} S^2, \qquad
   W\cdot W = \frac{1}{2} W^2,
   \quad
   S^\star = \sqrt{\max(S\cdot S, 0)}.
   \end{equation}

4. Compute the nonlinear-compatible $C_\mu$
   \begin{equation}
   C_\mu = Cmu_\text{nonlinear}(C_{a0},\dots, C_{a3}, \text{inv}, k, \epsilon),
   \end{equation}
   and then form the denominator
   \begin{equation}
   \text{denom} = (C_{NL6} + C_{NL7} S^\star) C_\mu.
   \end{equation}

5. From this denominator define
   \begin{equation}
   C_1 = \frac{C_{NL1}}{\text{denom}}, \quad
   C_2 = \frac{C_{NL2}}{\text{denom}}, \quad
   C_3 = \frac{C_{NL3}}{\text{denom}}.
   \end{equation}

6. Build the *identity tensor* $I$ and the quadratic tensor products:
   \begin{equation}
   SS = S S, \qquad WW = W W,
   \qquad WS = W S, \qquad SW^T = S W^T.
   \end{equation}

7. Construct the *quadratic bracket*:
   \begin{equation}
   \mathcal{Q} =
   C_1 \left( S S - \frac{S\cdot S}{3} I \right)
 + C_2 (W S + S W^T)
 + C_3 \left( W W - \frac{W\cdot W}{3} I \right).
   \end{equation}

8. The *quadratic nonlinear Reynolds-stress* contribution is
   \begin{equation}
   T_\text{quad} = -4 \mu_t \frac{k}{\epsilon} \, \mathcal{Q}.
   \end{equation}

   - If `model == Quadratic`, this is the final tensor returned.

9. For the *cubic* model, additional terms are included:

   - Compute
     \begin{equation}
     SSw = (S S) W,\qquad
     W^T S S = W^T (S S).
     \end{equation}
   - Define
     \begin{equation}
     \text{term}_4 = S S W + W^T S S,
     \end{equation}
     and
     \begin{equation}
     \text{term}_5 = (S\cdot S - W\cdot W) (S - W^T).
     \end{equation}

   - Let
     \begin{equation}
     C_4 = C_{NL4} C_\mu^2, \qquad
     C_5 = C_{NL5} C_\mu^2,
     \end{equation}
     and assemble the cubic bracket
     \begin{equation}
     \mathcal{C} = C_4 \, \text{term}_4 + C_5 \, \text{term}_5.
     \end{equation}

   - With $k^2/\epsilon^2 = k^2_\epsilon$ we obtain the cubic contribution
     \begin{equation}
     T_\text{cubic} = -8 \, \mu_t \, \frac{k^2}{\epsilon^2} \, \mathcal{C}.
     \end{equation}

10. In the cubic case, the total nonlinear stress tensor returned is $T_\text{NL} = T_\text{quad} + T_\text{cubic}$.


### 17.2 `computeGnl`

This function returns the *nonlinear production term* $G_\text{nl}$ associated with the
nonlinear Reynolds stress:

1. If `model == None`, it returns 0.
2. Otherwise it first computes $T_\text{NL}$ via `computeTRANS_NL`.
3. The production term is then
   \begin{equation}
   G_\text{nl} = T_\text{NL} : \nabla \mathbf{u},
   \end{equation}
   i.e. the double contraction between the nonlinear stress tensor and the velocity-gradient
   tensor.

This contribution is added to the standard production $G_k$ when nonlinear models are
enabled.


## 18. `computeCurvatureFactor`

This function computes a *curvature/rotation correction factor* $f_c$ (of Spalart–Shur
type) that multiplies the shear-production term.

1. If `model == None`, it returns $f_c = 1$.

2. Otherwise it starts from the invariants (clamped to non-negative):
   \begin{equation}
   S^2 = \max(\text{inv.S2}, 0), \qquad
   W^2 = \max(\text{inv.W2}, 0).
   \end{equation}

3. If $S^2 \le 0$, it again returns 1. Otherwise define
   \begin{equation}
   S\cdot S = \frac{1}{2} S^2, \qquad
   W\cdot W = \frac{1}{2} W^2.
   \end{equation}

4. Construct dimensionless measures of rotation vs strain:
   \begin{equation}
   r_* = \frac{W\cdot W}{\max(S\cdot S, 10^{-16})}, \qquad
   \tilde{r} = \frac{W\cdot W}{\max(S\cdot S + W\cdot W, 10^{-16})}.
   \end{equation}

5. With internal constants
   \begin{equation}
   C_\text{rot1} = 1.0,\quad C_\text{rot2} = 2.0,\quad
   C_\text{rot3} = 1.0,\quad C_\text{max} = 1.25,
   \end{equation}
   it forms a Spalart–Shur-like correction
   \begin{equation}
   f_\text{rot} =
   (1 + C_\text{rot1}) \,
   \frac{2 \tilde{r}}{1 + \tilde{r}}
   \left( 1 - C_\text{rot2} \frac{r_*}{1 + C_\text{rot3} r_*} \right)
   - C_\text{rot1}.
   \end{equation}

6. The function then clamps $f_\text{rot}$ to be at least $-1$ and finally defines
   \begin{equation}
   f_c = \min(C_\text{max}, 1 + f_\text{rot}).
   \end{equation}

The resulting $f_c$ is used to scale production in curvature/rotation-corrected
$k\text{–}\epsilon$ models, enhancing or damping turbulence depending on the local ratio
of rotation to strain.


## Corrections implemented in the $k$–$\epsilon$ models

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
