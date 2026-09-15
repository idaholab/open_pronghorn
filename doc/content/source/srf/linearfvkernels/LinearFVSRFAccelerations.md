# LinearFVSRFAccelerations

!syntax description /LinearFVKernels/LinearFVSRFAccelerations

## Description

`LinearFVSRFAccelerations` adds the apparent acceleration terms that arise when the
momentum equations are solved in a rotating body reference frame. For momentum component
$i$, the elemental right-hand-side contribution is

\begin{equation}
b_i = -\rho \left[
\boldsymbol{\omega}_B \times
  \left(\boldsymbol{\omega}_B \times \boldsymbol{r}_{mc}\right)
+ \dot{\boldsymbol{\omega}}_B \times \boldsymbol{r}_{mc}
+ 2\boldsymbol{\omega}_B \times \boldsymbol{u}_B
\right]_i.
\end{equation}

The three terms inside the brackets are the centripetal, Euler, and Coriolis acceleration
terms, respectively. The leading minus sign supplies the corresponding apparent forces
in the rotating-frame momentum equation.

The density, angular velocity, angular acceleration, and metacenter displacement are
provided as functors. [`LinearFVSRFFunctorMaterial.md`] creates the SRF kinematic functors
expected by this kernel. The `u`, `v`, and optional `w` parameters identify the body-frame
velocity components used in the Coriolis term. If `w` is omitted, its value is taken as
zero.

One kernel is required for each solved momentum component. As with other linear FV body
forces, include the names of these kernels in the `body_force` parameter of the
Rhie-Chow interpolator when pressure-velocity coupling is used.

## Example input syntax

The following example adds the SRF acceleration source to the $x$-momentum equation:

```text
[LinearFVKernels]
  [u_srf_acceleration]
    type = LinearFVSRFAccelerations
    variable = vel_x
    momentum_component = x
    rho = rho
    u = vel_x
    v = vel_y
    w = vel_z
    omega_brf = omega_brf
    omega_dot_brf = omega_dot_brf
    r_mc = r_mc
  []
[]
```

Create analogous objects with `momentum_component = y` and, for a three-dimensional
model, `momentum_component = z`.

!syntax parameters /LinearFVKernels/LinearFVSRFAccelerations

!syntax inputs /LinearFVKernels/LinearFVSRFAccelerations
