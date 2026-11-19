# LinearFVTurbulentBlendedDiffusion

`LinearFVTurbulentBlendedDiffusion` is the diffusion kernel [LinearFVTurbulentDiffusion.md] with
each matrix row and right-hand-side entry multiplied by a per-element weight. The weight is 0 in
wall-bounded elements, the
[!param](/LinearFVKernels/LinearFVTurbulentBlendedDiffusion/blending_factor) $\lambda$ elsewhere if
one is given, and 1 otherwise. The face discretization, including the optional non-orthogonal
correction, is unchanged.

See [LinearFVTurbulentBlendedAdvection.md] for:

- how the weight is defined and applied on internal and boundary faces,
- why every transport kernel of $\epsilon$ must use this kernel family in the two-layer
  k–$\epsilon$ variants.

Both diffusion terms of the $\epsilon$ equation, the molecular term and the turbulent
($\mu_t/\sigma_\epsilon$) term, must use this kernel with the same `blending_factor`:

!listing test/tests/kEpsilon/special-cases/channel_ERCOFTAC_buoyant.i block=LinearFVKernels/TKED_diffusion

!listing test/tests/kEpsilon/special-cases/channel_ERCOFTAC_buoyant.i block=LinearFVKernels/TKED_turb_diffusion

Without `blending_factor`, the kernel gives the same result as [LinearFVTurbulentDiffusion.md].

!syntax parameters /LinearFVKernels/LinearFVTurbulentBlendedDiffusion

!syntax inputs /LinearFVKernels/LinearFVTurbulentBlendedDiffusion

!syntax children /LinearFVKernels/LinearFVTurbulentBlendedDiffusion
