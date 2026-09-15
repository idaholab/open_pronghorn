# LinearFVSRFFunctorMaterial

!syntax description /FunctorMaterials/LinearFVSRFFunctorMaterial

## Description

`LinearFVSRFFunctorMaterial` defines the kinematic functors used by the single rotating
frame (SRF) linear finite-volume objects. It describes the position and motion of the
body reference frame relative to the metacenter reference frame.

The material provides two input modes:

- `fixed` uses prescribed pitch, yaw, and roll angles together with prescribed angular
  velocity and angular acceleration components.
- `pitch_yaw_roll` describes each angle with a sinusoid. For an angle $q$, the motion is

  \begin{equation}
  q(t) = A_q \sin\left(\frac{2\pi t}{T_q} + \delta_q\right),
  \end{equation}

  where $A_q$, $T_q$, and $\delta_q$ are its amplitude, period, and phase. Amplitudes
  and phases are supplied in degrees, periods are supplied in seconds, and the generated
  angle functors are in radians. Each period must be greater than zero.

The three sinusoidal motions are converted to body-frame angular velocity using

\begin{equation}
\boldsymbol{\omega}_B =
\begin{bmatrix}
\dot{\phi} + \sin(\psi)\dot{\theta} \\
\cos(\phi)\cos(\psi)\dot{\theta} + \sin(\phi)\dot{\psi} \\
-\sin(\phi)\cos(\psi)\dot{\theta} + \cos(\phi)\dot{\psi}
\end{bmatrix},
\end{equation}

where $\phi$, $\theta$, and $\psi$ denote roll, pitch, and yaw, respectively. Consequently,
when more than one angle varies, the components of `omega_brf` are not simply the three
Euler-angle rates. The material also differentiates this expression to construct
`omega_dot_brf`.

The displacement functor is

\begin{equation}
\boldsymbol{r}_{mc} = \boldsymbol{x} - \boldsymbol{x}_{mc},
\end{equation}

where $\boldsymbol{x}_{mc}$ is set with `mc_origin`.

## Generated functors

The following functors are available to other objects:

| Functor | Meaning | Units |
| :- | :- | :- |
| `r_mc` | Displacement from the metacenter origin to the evaluation point | length |
| `pitch_angle`, `yaw_angle`, `roll_angle` | Orientation angles | rad |
| `omega_brf` | Angular velocity in the body reference frame | rad/s |
| `omega_dot_brf` | Angular acceleration in the body reference frame | rad/s$^2$ |
| `omega_pitch`, `omega_yaw`, `omega_roll` | Scalar angular-velocity outputs | rad/s |
| `omega_dot_pitch`, `omega_dot_yaw`, `omega_dot_roll` | Scalar angular-acceleration outputs | rad/s$^2$ |

The vector functors can be passed directly to
[`LinearFVSRFAccelerations.md`], while the angle functors are consumed by
[`LinearFVSRFMomentumBoussinesq.md`] and [`LinearFVSRFSource.md`].

## Example input syntax

The following example defines a fixed orientation:

!listing test/tests/ocean_MRF/diff_heated_static_tilt.i block=FunctorMaterials/SRF_Functor_Material

!syntax parameters /FunctorMaterials/LinearFVSRFFunctorMaterial

!syntax inputs /FunctorMaterials/LinearFVSRFFunctorMaterial
