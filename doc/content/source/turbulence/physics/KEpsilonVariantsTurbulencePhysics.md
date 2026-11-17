# KEpsilonVariantsTurbulencePhysics

`KEpsilonVariantsTurbulencePhysics` sets up a k-$\epsilon$ turbulence model with a linear finite
volume discretization, for use alongside a segregated flow Physics such as
[WCNSLinearFVFlowPhysics.md]. It derives from MOOSE's
[WCNSLinearFVTurbulencePhysics.md] and reuses its variables, advection and diffusion kernels,
boundary conditions and effective viscosity materials.

With respect to the MOOSE Physics, which only implements the standard high-Reynolds-number
k-$\epsilon$ model, this Physics adds the model variants and corrections implemented in
OpenPronghorn. It does so by creating the OpenPronghorn [kEpsilonTKESourceSink.md],
[kEpsilonTKEDSourceSink.md] and [kEpsilonViscosity.md] objects in place of the MOOSE ones.

## Model variants

The variant is selected with the `k_epsilon_variant` parameter:

- `Standard`: high-Reynolds-number standard k-$\epsilon$,
- `StandardLowRe`: low-Reynolds-number standard k-$\epsilon$,
- `StandardTwoLayer`: two-layer standard k-$\epsilon$,
- `Realizable`: realizable k-$\epsilon$, with a variable $C_\mu$,
- `RealizableTwoLayer`: two-layer realizable k-$\epsilon$.

For the two-layer variants, `two_layer_flavor` selects the near-wall formulation
(`Wolfstein`, `NorrisReynolds` or `Xu`). See the [theory](theory/turbulenceModeling.md) section
for the formulation of each variant.

The following optional corrections can be enabled independently of the variant: `use_buoyancy`,
`use_compressibility`, `use_yap`, `use_low_re_Gprime`, `nonlinear_model` and `curvature_model`.
The underlying kernels error out when a correction is not applicable to the selected variant.

## Wall distance

The two-layer and low-Reynolds-number formulations, as well as the Yap and low-Reynolds extra
production corrections, need the distance to the nearest wall. When one of these is selected,
this Physics creates an auxiliary variable and a [WallDistanceAux.md] computing the distance to
the boundaries listed in `turbulence_walls`. To use an existing wall distance field instead, pass
its name with the `wall_distance` parameter, in which case no object is created.

## Coupling with other Physics

The buoyancy production term requires a temperature field, a thermal expansion coefficient and a
gravity vector. These are not specified on this Physics: the temperature is taken from the coupled
fluid heat transfer Physics, while `thermal_expansion` and `gravity` are taken from the coupled
flow Physics. Setting `use_buoyancy = true` without a fluid heat transfer Physics is an error.

!alert note
The OpenPronghorn k-$\epsilon$ objects take the `C1_eps`, `C2_eps` and `Pr_t` closure
coefficients as constants, while the base Physics declares them as functors. Specifying a
non-constant functor for one of these is an error.

!syntax parameters /Physics/NavierStokes/AdvancedKEpsilon/KEpsilonVariantsTurbulencePhysics

!syntax inputs /Physics/NavierStokes/AdvancedKEpsilon/KEpsilonVariantsTurbulencePhysics

!syntax children /Physics/NavierStokes/AdvancedKEpsilon/KEpsilonVariantsTurbulencePhysics
