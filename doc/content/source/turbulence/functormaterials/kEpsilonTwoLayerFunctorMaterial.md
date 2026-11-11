# kEpsilonTwoLayerFunctorMaterial

`kEpsilonTwoLayerFunctorMaterial` is a functor material that computes the near-wall quantities of
the two-layer k–$\epsilon$ treatment. Every object that uses the two-layer model reads them from
this one material: [kEpsilonViscosity.md], [kEpsilonTKEDSourceSink.md] and the blended transport
kernels ([LinearFVTurbulentBlendedAdvection.md], [LinearFVTurbulentBlendedDiffusion.md],
[LinearFVTurbulentBlendedTimeDerivative.md]). This guarantees that $\mu_t$ and the
$\epsilon$-equation use the same blending function and the same near-wall length scales.

!alert note
The explanations here are a summary. See the [theory](theory/turbulenceModeling.md) for more
details on the two-layer models.

## Computed functors

The material defines three functors. Their names are set with the
[!param](/FunctorMaterials/kEpsilonTwoLayerFunctorMaterial/two_layer_lambda),
[!param](/FunctorMaterials/kEpsilonTwoLayerFunctorMaterial/two_layer_epsilon) and
[!param](/FunctorMaterials/kEpsilonTwoLayerFunctorMaterial/two_layer_mu_t) parameters. They default
to `two_layer_lambda`, `two_layer_epsilon` and `two_layer_mu_t`.

All three are based on the wall-distance Reynolds number

\begin{equation}
Re_d = \frac{\sqrt{k}\, d}{\nu}, \qquad \nu = \frac{\mu}{\rho},
\end{equation}

where $d$ is the [!param](/FunctorMaterials/kEpsilonTwoLayerFunctorMaterial/wall_distance) functor,
typically computed with [WallDistanceAux.md]. Negative values of $k$ are clipped to zero.

### Blending factor

\begin{equation}
\lambda = \tfrac12\left[1 + \tanh\!\left(\frac{Re_d - Re^*}{A}\right)\right],
\qquad A = \frac{\Delta Re}{\operatorname{atanh}(0.98)},
\end{equation}

with $Re^*$ = [!param](/FunctorMaterials/kEpsilonTwoLayerFunctorMaterial/Rey_star) and
$\Delta Re$ = [!param](/FunctorMaterials/kEpsilonTwoLayerFunctorMaterial/DeltaRey). This gives
$\lambda \to 0$ in the viscosity-affected region ($Re_d \ll Re^*$) and $\lambda \to 1$ in the
fully turbulent region. $\lambda$ goes from 0.01 to 0.99 over $Re^* \pm \Delta Re$.

### Algebraic near-wall dissipation

\begin{equation}
\epsilon_{2L} = \frac{k^{3/2}}{\ell_\epsilon}.
\end{equation}

### Algebraic near-wall turbulent viscosity

\begin{equation}
\mu_{t,2L} = \left(\frac{\mu_t}{\mu}\right)_{2L} \mu.
\end{equation}

## Two-layer flavors

The dissipation length $\ell_\epsilon$ and the viscosity ratio $(\mu_t/\mu)_{2L}$ depend on
[!param](/FunctorMaterials/kEpsilonTwoLayerFunctorMaterial/two_layer_flavor). With
$c_\ell = 0.42\, C_\mu^{-3/4}$:

- `Wolfstein` [!cite](wolfstein1969velocity):
  \begin{equation}
  \ell_\epsilon = c_\ell\, d \left(1 - e^{-Re_d/(2 c_\ell)}\right), \qquad
  \left(\frac{\mu_t}{\mu}\right)_{2L} = 0.42\, C_\mu^{1/4} Re_d \left(1 - e^{-Re_d/70}\right).
  \end{equation}
- `NorrisReynolds` [!cite](norris1975oneequation):
  \begin{equation}
  \ell_\epsilon = c_\ell\, d \frac{Re_d}{Re_d + 5.3}, \qquad
  \left(\frac{\mu_t}{\mu}\right)_{2L} = 0.42\, C_\mu^{1/4} Re_d \left(1 - e^{-Re_d/50.5}\right).
  \end{equation}

For both flavors,
$(\mu_t/\mu)_{2L} = \rho C_\mu \ell_\mu \sqrt{k}/\mu$ with $\ell_\mu = c_\ell\, d\, (1 - e^{-Re_d/A_\mu})$.
This is the one-equation viscosity of the Chen–Patel two-layer model.

## Usage

Declare the material once, and point the two-layer parameters of the other objects at its outputs.
The kernels check that the functors they receive are outputs of a
`kEpsilonTwoLayerFunctorMaterial`. Every transport kernel of the $\epsilon$ equation must also be a
blended kernel that uses the same `two_layer_lambda`.

!listing test/tests/kEpsilon/special-cases/channel_ERCOFTAC_buoyant.i block=FunctorMaterials/two_layer

!listing test/tests/kEpsilon/special-cases/channel_ERCOFTAC_buoyant.i block=LinearFVKernels/TKED_source_sink

!syntax parameters /FunctorMaterials/kEpsilonTwoLayerFunctorMaterial

!syntax inputs /FunctorMaterials/kEpsilonTwoLayerFunctorMaterial

!syntax children /FunctorMaterials/kEpsilonTwoLayerFunctorMaterial
