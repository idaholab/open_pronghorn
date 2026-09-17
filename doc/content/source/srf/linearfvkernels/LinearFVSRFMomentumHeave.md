# LinearFVSRFMomentumHeave

!syntax description /LinearFVKernels/LinearFVSRFMomentumHeave

## Description

`LinearFVSRFMomentumHeave` adds the apparent acceleration associated with harmonic heave
motion to a linear finite-volume momentum equation. Unlike a Boussinesq source, the
heave source depends only on density and the prescribed frame acceleration. It does not
contain a temperature difference or a thermal expansion coefficient.

The prescribed heave displacement is

\begin{equation}
z_h = A_h \sin\left(\omega_h t + \delta_h\right),
\end{equation}

where $A_h$ is `amplitude`, $\delta_h$ is `phase`, and

\begin{equation}
\omega_h = \frac{2\pi}{T_h},
\end{equation}

with $T_h$ given by `period`. The acceleration of the moving frame is

\begin{equation}
\ddot{z}_h = -A_h\omega_h^2
\sin\left(\omega_h t + \delta_h\right).
\end{equation}

The apparent acceleration acting on the fluid is opposite to the frame acceleration:

\begin{equation}
a_h = -\ddot{z}_h = A_h\omega_h^2
\sin\left(\omega_h t + \delta_h\right).
\end{equation}

Consequently, the elemental right-hand-side source density for the momentum component
selected by `momentum_component` is

\begin{equation}
b_i = \rho a_h.
\end{equation}

For the standard ship coordinate convention, heave acts in the vertical $z$ direction,
so the kernel is normally assigned to the $z$-momentum equation with
`momentum_component = z`. The kernel applies the heave acceleration directly to the
selected body-frame component; it does not rotate the translational acceleration using
pitch, yaw, or roll angles.

The density is supplied as a functor and can therefore vary in space and time. The
period must be greater than zero. As with other linear FV body forces, include the name
of this kernel in the `body_force` parameter of the Rhie-Chow interpolator when
pressure-velocity coupling is used.

## Example input syntax

The following example applies harmonic heave to the vertical momentum equation:

```text
[LinearFVKernels]
  [z_heave]
    type = LinearFVSRFMomentumHeave
    variable = vel_z
    momentum_component = z
    rho = rho
    amplitude = 5.7
    period = 9.4
    phase = 0
  []
[]
```

!syntax parameters /LinearFVKernels/LinearFVSRFMomentumHeave

!syntax inputs /LinearFVKernels/LinearFVSRFMomentumHeave
