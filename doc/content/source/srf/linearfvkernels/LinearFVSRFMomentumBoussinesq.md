# LinearFVSRFMomentumBoussinesq

!syntax description /LinearFVKernels/LinearFVSRFMomentumBoussinesq

## Description

`LinearFVSRFMomentumBoussinesq` supplies a Boussinesq buoyancy source to a momentum
equation solved in the body reference frame. The user specifies gravity in the inertial
metacenter frame. The kernel rotates it into the body frame before evaluating the source:

\begin{equation}
\boldsymbol{g}_B = \boldsymbol{R}_{I\rightarrow B}
  (\theta,\psi,\phi)\boldsymbol{g}_I,
\end{equation}

\begin{equation}
b_i = -\rho\alpha\left(T-T_{ref}\right)(\boldsymbol{g}_B)_i,
\end{equation}

where $\theta$, $\psi$, and $\phi$ are pitch, yaw, and roll; $\rho$ is the reference
density; $\alpha$ is the thermal expansion coefficient; and $T_{ref}$ is the reference
temperature. The rotation follows the convention documented in [`SRFUtils.md`].

The angle inputs are functors and must evaluate in radians. They can vary in space and
time, although the usual SRF setup obtains them from
[`LinearFVSRFFunctorMaterial.md`]. One kernel is required for each solved momentum
component.

## Example input syntax

This example applies the rotated Boussinesq source to the $x$-momentum equation:

!listing test/tests/ocean_MRF/diff_heated_static_tilt.i block=LinearFVKernels/u_buoyancy

!syntax parameters /LinearFVKernels/LinearFVSRFMomentumBoussinesq

!syntax inputs /LinearFVKernels/LinearFVSRFMomentumBoussinesq
