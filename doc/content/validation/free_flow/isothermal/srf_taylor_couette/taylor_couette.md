# Taylor--Couette Flow in a Single Rotating Frame

!tag name=Rotating Frame of Reference: Taylor Couette
    image=../../media/validation/free_flow/isothermal/ercoftac_030_bfs/titlepage_srftaylorcouette.pdf
    description=Non-Inertial Frame vs Inertial Frame formulation of a Taylor Couette flow in a cylinder
    pairs=compressibility:incompressible
          heattransfer:isothermal
          convection_type:forced
          transient:steady
          flow_regime:laminar
          fluid:air
          flow_configuration:free-flow
          number_of_phases:one

## Problem Description

This validation problem assesses the single-rotating-frame (SRF) formulation in
`OpenPronghorn` using steady, laminar Taylor--Couette flow between two concentric
cylinders. The detailed description of the physics can be found in [!cite](taylor1923stability).
In the inertial frame, the inner cylinder rotates at a prescribed
angular velocity while the outer cylinder remains stationary. The computational
frame rotates with the inner cylinder so that the inner wall is stationary in
the rotating frame and the outer wall moves with the corresponding relative
angular velocity.

The inertial-frame configuration is defined by an inner radius $R_i=0.35$ m,
an outer radius $R_o=1.0$ m, an inner-cylinder angular velocity
$\Omega_i=-0.01$ rad/s, and a stationary outer cylinder,
$\Omega_o=0$. The computational-frame angular velocity is
$\Omega_f=\Omega_i$. Therefore,

!equation id=tc-relative-omega
\Omega_{i,\mathrm{rel}} = \Omega_i-\Omega_f = 0,
\qquad
\Omega_{o,\mathrm{rel}} = \Omega_o-\Omega_f = 0.01\ \mathrm{rad/s}.

The purpose of this case is to verify that the velocity and pressure fields
computed in the rotating frame reproduce the analytical inertial-frame
Taylor--Couette solution after the velocity is transformed back to the inertial
frame.

## Modeling Parameters

The fluid is incompressible and Newtonian. The parameters used for the
validation case are summarized below.

| Parameter | Value | Units |
| :- | :- | :- |
| Inner-cylinder radius, $R_i$ | 0.35 | m |
| Outer-cylinder radius, $R_o$ | 1.0 | m |
| Density, $\rho$ | 1.0 | kg/m$^3$ |
| Dynamic viscosity, $\mu$ | $5.0\times10^{-4}$ | Pa s |
| Inner-cylinder angular velocity, $\Omega_i$ | -0.01 | rad/s |
| Outer-cylinder angular velocity, $\Omega_o$ | 0.0 | rad/s |
| Rotating-frame angular velocity, $\Omega_f$ | -0.01 | rad/s |
| Pressure-reference radius | 0.75 | m |

For this low-Reynolds-number configuration, the steady laminar solution is
axisymmetric and has no radial velocity. The analytical tangential velocity is

!equation id=tc-velocity
u_\theta(r) = A r + \frac{B}{r},

where

!equation id=tc-coeff-a
A = \frac{\Omega_o R_o^2-\Omega_i R_i^2}{R_o^2-R_i^2},

and

!equation id=tc-coeff-b
B = \frac{R_i^2R_o^2\left(\Omega_i-\Omega_o\right)}{R_o^2-R_i^2}.

The corresponding radial pressure gradient follows from radial equilibrium,

!equation id=tc-pressure-gradient
\frac{dp}{dr} = \rho\frac{u_\theta^2}{r}.

Using [!eqref](tc-velocity), the analytical pressure distribution can be written
apart from an arbitrary additive constant as

!equation id=tc-pressure
p(r) = \rho\left[
\frac{A^2}{2}r^2
+2AB\ln(r)
-\frac{B^2}{2r^2}
\right] + C.

Because incompressible pressure is defined only up to an additive constant,
the numerical and analytical pressure profiles are shifted to zero at the
sampled finite-volume centroid closest to $r=0.75$ m before comparison.

## `OpenPronghorn` Model

The annular domain is discretized with a two-dimensional finite-volume mesh.
The solution is obtained with the linear SIMPLE solver. The momentum equations
include the SRF acceleration terms through `LinearFVSRFAccelerations`, while
`LinearFVSRFFunctorMaterial` provides the prescribed frame angular velocity and
rotation-center quantities. The inner wall is stationary in the rotating frame
and the outer wall is assigned the relative tangential wall velocity.

The velocity solved by the model is the rotating-frame velocity,
$\mathbf{u}_{\mathrm{rel}}$. For comparison with the analytical inertial-frame
solution, it is reconstructed according to

!equation id=tc-velocity-transform
\mathbf{u}_{\mathrm{abs}}
=
\mathbf{u}_{\mathrm{rel}}
+
\boldsymbol{\Omega}_f\times\mathbf{r}.

For rotation about the $z$ axis, the reconstructed Cartesian components are

!equation id=tc-ux-transform
u_{\mathrm{abs},x}
=
u_{\mathrm{rel},x}-\Omega_f y,

and

!equation id=tc-uy-transform
u_{\mathrm{abs},y}
=
u_{\mathrm{rel},y}+\Omega_f x.

An `ElementValueSampler` records `vel_abs_x`, `vel_abs_y`, and `pressure` at the
finite-volume element centroids. The validation script selects the radial row
of centroids closest to the positive $x$ axis and evaluates the analytical
solution at those exact centroid radii. For the current mesh, the selected row
is located at approximately $0.15^\circ$ from the positive $x$ axis. Along this
row, `vel_abs_y` is used as the tangential velocity and `vel_abs_x` as the
radial velocity; the angular offset is sufficiently small that its effect on
the tangential-velocity comparison is negligible.

The input file for the SRF Taylor--Couette calculation is embedded below.

!listing validation/free_flow/isothermal/srf_taylor_couette/taylor_couette_2d_rel.i

## Results

The reconstructed inertial-frame tangential velocity is compared directly with
[!eqref](tc-velocity) at the same finite-volume centroid radii. The numerical
profile closely follows the analytical Taylor--Couette solution throughout the
annular gap.

The validation figures are generated during the MooseDocs build using
`taylor_couette_plot.py`. The plotting script reads the accepted numerical
radial-profile data from `taylor_couette_results.csv` and evaluates the
analytical velocity and pressure solutions at the same radial locations.

!media media/validation/free_flow/isothermal/srf_taylor_couette/taylor_couette_plot.py
       image_name=taylor_couette_velocity_comparison.png id=fig:taylor-couette-velocity
       style=width:80%;margin-left:auto;margin-right:auto;text-align:center
       caption=Reconstructed inertial-frame tangential velocity compared with the analytical Taylor--Couette solution.

The numerical pressure is also compared with the analytical pressure profile
from [!eqref](tc-pressure). Both pressure profiles are shifted using the same
reference centroid because only pressure differences are physically relevant
for this incompressible problem.

!media media/validation/free_flow/isothermal/srf_taylor_couette/taylor_couette_plot.py
       image_name=taylor_couette_pressure_comparison.png id=fig:taylor-couette-pressure
       style=width:80%;margin-left:auto;margin-right:auto;text-align:center
       caption=Shifted numerical pressure profile compared with the analytical Taylor--Couette pressure distribution.

For the current validation calculation, the reconstructed tangential velocity
has an RMS normalized error of approximately $0.29\%$ and a maximum normalized
error of approximately $0.47\%$. The shifted pressure profile has an RMS error
of approximately $0.18\%$ of the analytical pressure span and a maximum error
of approximately $0.52\%$ of the analytical pressure span.

## Validation

The velocity validation metric is the pointwise difference between the
reconstructed inertial-frame tangential velocity and the analytical solution,
normalized by the physical inner-wall speed,

!equation id=tc-velocity-error
\epsilon_u(r_i)
=
\frac{
\left|u_{\theta,\mathrm{OP}}(r_i)-u_{\theta,\mathrm{analytic}}(r_i)\right|
}{
\left|\Omega_i R_i\right|
}.

The pressure validation metric is the pointwise shifted-pressure error
normalized by the total analytical pressure variation over the sampled radial
profile,

!equation id=tc-pressure-error
\epsilon_p(r_i)
=
\frac{
\left|p_{\mathrm{OP}}'(r_i)-p_{\mathrm{analytic}}'(r_i)\right|
}{
\max\left(p_{\mathrm{analytic}}'\right)-
\min\left(p_{\mathrm{analytic}}'\right)
},

where the prime denotes pressure shifted by its value at the common reference
centroid.

The automated validation requires both pointwise error measures to remain below
$1\%$:

!equation id=tc-validation-bounds
\epsilon_u \le 0.01,
\qquad
\epsilon_p \le 0.01.

The current solution satisfies both validation criteria. The maximum velocity
error is approximately $0.47\%$ and the maximum pressure error is approximately
$0.52\%$ using the normalizations defined above. The validation script also
writes `taylor_couette_analytical_comparison.csv`, which contains the sampled
centroid coordinates, numerical and analytical profiles, and the associated
pointwise errors.

The validation implementation is embedded below.

!listing validation/free_flow/isothermal/srf_taylor_couette/taylor_couette_vnv.py
         id=tc-vnv
         caption=Validation script for the SRF Taylor--Couette velocity and pressure profiles.
