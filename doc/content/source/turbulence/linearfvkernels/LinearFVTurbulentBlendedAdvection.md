# LinearFVTurbulentBlendedAdvection

`LinearFVTurbulentBlendedAdvection` is the advection kernel [LinearFVTurbulentAdvection.md] with
each matrix row and right-hand-side entry multiplied by a per-element weight. The face
discretization itself is unchanged.

It belongs to a family of blended turbulence transport kernels:

- `LinearFVTurbulentBlendedAdvection` (this kernel), based on [LinearFVTurbulentAdvection.md],
- [LinearFVTurbulentBlendedDiffusion.md], based on [LinearFVTurbulentDiffusion.md],
- [LinearFVTurbulentBlendedTimeDerivative.md], based on [LinearFVTurbulentTimeDerivative.md].

## Row weight

For the row of element $e$, the weight is

\begin{equation}
w_e =
\begin{cases}
0 & \text{if } e \text{ is wall bounded (touches a boundary in } \texttt{walls}\text{)}, \\
\lambda_e & \text{if a } \texttt{blending\_factor} \text{ is provided}, \\
1 & \text{otherwise.}
\end{cases}
\end{equation}

Without [!param](/LinearFVKernels/LinearFVTurbulentBlendedAdvection/blending_factor), the kernel
gives the same result as [LinearFVTurbulentAdvection.md], which skips wall-bounded rows.

## Use with the two-layer k–$\epsilon$ model

For the `StandardTwoLayer` and `RealizableTwoLayer` variants, [kEpsilonTKEDSourceSink.md] blends
the $\epsilon$-row between the transport equation and the algebraic near-wall value
$\epsilon_{2L}$:

\begin{equation}
\lambda \left(\text{transport} - \text{sources}\right) + (1 - \lambda)\, a_{2L}
\left(\epsilon - \epsilon_{2L}\right) = 0,
\qquad a_{2L} = \frac{\rho\, \epsilon_{2L}}{k}.
\end{equation}

The source/sink kernel only provides its own part of that row. Every transport kernel acting on
$\epsilon$ (advection, diffusion and, for transient runs, the time derivative) must therefore be a
blended kernel. Its `blending_factor` must be the `two_layer_lambda` output of the
[kEpsilonTwoLayerFunctorMaterial.md]. [kEpsilonTKEDSourceSink.md] checks this at setup. It raises
an error if any transport kernel on $\epsilon$ is not blended, or uses a different blending factor.

!listing test/tests/kEpsilon/special-cases/channel_ERCOFTAC_buoyant.i block=LinearFVKernels/TKED_advection

!alert note
Use the blended kernels only for the $\epsilon$ equation. The $k$ equation is solved in all cells,
so its transport kernels keep the regular [LinearFVTurbulentAdvection.md] and
[LinearFVTurbulentDiffusion.md] kernels.

!syntax parameters /LinearFVKernels/LinearFVTurbulentBlendedAdvection

!syntax inputs /LinearFVKernels/LinearFVTurbulentBlendedAdvection

!syntax children /LinearFVKernels/LinearFVTurbulentBlendedAdvection
