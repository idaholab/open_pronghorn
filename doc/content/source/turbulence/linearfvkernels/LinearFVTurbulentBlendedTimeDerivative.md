# LinearFVTurbulentBlendedTimeDerivative

`LinearFVTurbulentBlendedTimeDerivative` is the time-derivative kernel
[LinearFVTurbulentTimeDerivative.md] with its diagonal and right-hand-side entries multiplied by
a per-element weight. The weight is:

- 0 in wall-bounded elements,
- the [!param](/LinearFVKernels/LinearFVTurbulentBlendedTimeDerivative/blending_factor) $\lambda$
  elsewhere, if one is given,
- 1 otherwise.

In transient runs with the `StandardTwoLayer` or `RealizableTwoLayer` k–$\epsilon$ variants, use
this kernel for the $\epsilon$ equation with `blending_factor` set to the `two_layer_lambda` output
of [kEpsilonTwoLayerFunctorMaterial.md]. Then, in the viscosity-affected region ($\lambda \to 0$),
the time term does not pull $\epsilon$ back toward its old value, and the row relaxes toward the
algebraic $\epsilon_{2L}$. See [LinearFVTurbulentBlendedAdvection.md] for how the blended row is
built and for the consistency check done by [kEpsilonTKEDSourceSink.md].

Without `blending_factor`, the kernel gives the same result as
[LinearFVTurbulentTimeDerivative.md].

!syntax parameters /LinearFVKernels/LinearFVTurbulentBlendedTimeDerivative

!syntax inputs /LinearFVKernels/LinearFVTurbulentBlendedTimeDerivative

!syntax children /LinearFVKernels/LinearFVTurbulentBlendedTimeDerivative
