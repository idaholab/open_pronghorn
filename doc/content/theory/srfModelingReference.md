# Single Rotating Frame Equation and Implementation Reference

This page is a compact reference for the non-inertial-frame equations implemented by the
OpenPronghorn SRF objects. The full derivation and physical interpretation are provided
in [`srfModeling.md`].

## Coordinate convention

| Quantity | Definition |
| :- | :- |
| Frame $I$ | Inertial frame in which gravity and environmental vectors are prescribed |
| Frame $B$ | Vessel- and reactor-fixed computational frame |
| $x_B$ | Longitudinal axis; roll axis |
| $y_B$ | Transverse axis; pitch axis |
| $z_B$ | Vertical axis; yaw axis and conventional heave direction |
| $\phi$ | Roll angle |
| $\theta$ | Pitch angle |
| $\psi$ | Yaw angle |
| $\mathbf r_B$ | Vector from `mc_origin` to the cell center, expressed in frame $B$ |
| $\mathbf u_B$ | Fluid velocity relative to the body-fixed mesh |

The transformation directions are

\begin{equation}
\mathbf v_I=\mathbf C_{IB}\mathbf v_B,
\qquad
\mathbf v_B=\mathbf C_{BI}\mathbf v_I,
\qquad
\mathbf C_{BI}=\mathbf C_{IB}^{T}.
\end{equation}

## Rotation matrices

The implemented body-to-inertial sequence is

\begin{equation}
\mathbf C_{IB}
=\mathbf R_y\theta\,\mathbf R_z\psi\,\mathbf R_x\phi.
\end{equation}

The inertial-to-body matrix is

\begin{equation}
\mathbf C_{BI}=
\begin{bmatrix}
c_\psi c_\theta & s_\psi & -c_\psi s_\theta \\
s_\phi s_\theta-c_\phi s_\psi c_\theta & c_\phi c_\psi &
s_\phi c_\theta+c_\phi s_\psi s_\theta \\
c_\phi s_\theta+s_\phi s_\psi c_\theta & -s_\phi c_\psi &
c_\phi c_\theta-s_\phi s_\psi s_\theta
\end{bmatrix},
\end{equation}

where $c_q=\cos q$ and $s_q=\sin q$.

## Angular kinematics

For the implemented pitch-yaw-roll convention,

\begin{equation}
\boldsymbol{\omega}_B=
\begin{bmatrix}
\dot\phi+s_\psi\dot\theta \\
c_\phi c_\psi\dot\theta+s_\phi\dot\psi \\
-s_\phi c_\psi\dot\theta+c_\phi\dot\psi
\end{bmatrix}.
\end{equation}

The angular acceleration is

\begin{equation}
\dot{\boldsymbol{\omega}}_B=
\begin{bmatrix}
\ddot\phi+s_\psi\ddot\theta+c_\psi\dot\psi\dot\theta \\
c_\phi c_\psi\ddot\theta+s_\phi\ddot\psi
-s_\phi c_\psi\dot\phi\dot\theta
-c_\phi s_\psi\dot\psi\dot\theta
+c_\phi\dot\phi\dot\psi \\
-s_\phi c_\psi\ddot\theta+c_\phi\ddot\psi
-c_\phi c_\psi\dot\phi\dot\theta
+s_\phi s_\psi\dot\psi\dot\theta
-s_\phi\dot\phi\dot\psi
\end{bmatrix}.
\end{equation}

For each harmonic angle $q$,

\begin{equation}
q=A_q\sin\left(\Omega_qt+\delta_q\right),
\qquad
\dot q=A_q\Omega_q\cos\left(\Omega_qt+\delta_q\right),
\qquad
\ddot q=-A_q\Omega_q^2\sin\left(\Omega_qt+\delta_q\right),
\end{equation}

with $\Omega_q=2\pi/T_q$.

## Body-frame momentum equation

The body-frame momentum equation is

\begin{equation}
\frac{\partial\rho\mathbf u_B}{\partial t}
+\nabla_B\cdot\left(\rho\mathbf u_B\otimes\mathbf u_B\right)
=-\nabla_Bp+\nabla_B\cdot\boldsymbol\tau_B
+\rho\mathbf g_B+\mathbf s_B+\mathbf f_{NI,B},
\end{equation}

with

\begin{equation}
\mathbf f_{NI,B}
=-\rho\left[
\mathbf a_{O,B}
+\dot{\boldsymbol\omega}_B\times\mathbf r_B
+\boldsymbol\omega_B\times
  \left(\boldsymbol\omega_B\times\mathbf r_B\right)
+2\boldsymbol\omega_B\times\mathbf u_B
\right].
\end{equation}

## Source-term summary

| Contribution | Source density added to the body-frame momentum equation | OpenPronghorn object |
| :- | :- | :- |
| Centrifugal | $-\rho\boldsymbol\omega_B\times\left(\boldsymbol\omega_B\times\mathbf r_B\right)$ | [`LinearFVSRFAccelerations.md`] |
| Euler | $-\rho\dot{\boldsymbol\omega}_B\times\mathbf r_B$ | [`LinearFVSRFAccelerations.md`] |
| Coriolis | $-2\rho\boldsymbol\omega_B\times\mathbf u_B$ | [`LinearFVSRFAccelerations.md`] |
| Translational frame motion | $-\rho\mathbf a_{O,B}$ | Specialized as harmonic heave by [`LinearFVSRFMomentumHeave.md`] |
| Gravity transformation | $\mathbf g_B=\mathbf C_{BI}\mathbf g_I$ | [`SRFUtils.md`] |
| Boussinesq buoyancy | $-\rho_{ref}\alpha\left(T-T_{ref}\right)\mathbf g_B$ | [`LinearFVSRFMomentumBoussinesq.md`] |
| General vector source | $\gamma\mathbf C_{BI}\mathbf s_I$ | [`LinearFVSRFSource.md`] |

## Gravity components

For $\mathbf g_I=[0,0,-g]^T$,

\begin{equation}
\mathbf g_B=
\begin{bmatrix}
g c_\psi s_\theta \\
-g\left(s_\phi c_\theta+c_\phi s_\psi s_\theta\right) \\
-g\left(c_\phi c_\theta-s_\phi s_\psi s_\theta\right)
\end{bmatrix}.
\end{equation}

Useful limiting cases are

\begin{equation}
\theta=\psi=\phi=0:
\qquad
\mathbf g_B=
\begin{bmatrix}
0 \\ 0 \\ -g
\end{bmatrix},
\end{equation}

\begin{equation}
\psi=\phi=0:
\qquad
\mathbf g_B=
\begin{bmatrix}
g\sin\theta \\ 0 \\ -g\cos\theta
\end{bmatrix},
\end{equation}

and


\begin{equation}
\theta=\psi=0:
\qquad
\mathbf g_B=
\begin{bmatrix}
0 \\ -g\sin\phi \\ -g\cos\phi
\end{bmatrix}.
\end{equation}

Pure yaw leaves the gravity vector unchanged.

## Harmonic heave sign convention

The prescribed displacement, frame acceleration, and apparent acceleration are

\begin{equation}
z_h=A_h\sin\left(\Omega_ht+\delta_h\right),
\end{equation}

\begin{equation}
\ddot z_h=-A_h\Omega_h^2\sin\left(\Omega_ht+\delta_h\right),
\end{equation}

and

\begin{equation}
a_h=-\ddot z_h
=A_h\Omega_h^2\sin\left(\Omega_ht+\delta_h\right).
\end{equation}

The heave momentum source is

\begin{equation}
b_{h,i}=\rho a_h.
\end{equation}

There is no factor of $\alpha$, $T-T_{ref}$, or gravity in this source.

## Generated SRF functors

| Functor | Mathematical quantity | Units |
| :- | :- | :- |
| `r_mc` | $\mathbf r_B$ | length |
| `pitch_angle` | $\theta$ | rad |
| `yaw_angle` | $\psi$ | rad |
| `roll_angle` | $\phi$ | rad |
| `omega_brf` | $\boldsymbol\omega_B$ | rad/s |
| `omega_dot_brf` | $\dot{\boldsymbol\omega}_B$ | rad/s$^2$ |
| `omega_pitch`, `omega_yaw`, `omega_roll` | Scalar angular-velocity outputs | rad/s |
| `omega_dot_pitch`, `omega_dot_yaw`, `omega_dot_roll` | Scalar angular-acceleration outputs | rad/s$^2$ |

These functors are created by [`LinearFVSRFFunctorMaterial.md`].

## Diagnostic checks

Use the following checks when constructing or reviewing an SRF input:

1. $\mathbf C_{BI}=\mathbf I$ at zero pitch, yaw, and roll.
2. $\mathbf C_{IB}\mathbf C_{BI}=\mathbf I$ for every angle combination.
3. Pure yaw does not change $\mathbf g_B$.
4. The centrifugal and Euler terms vanish at `mc_origin`.
5. The Coriolis term vanishes when $\mathbf u_B=\mathbf 0$ or
   $\boldsymbol\omega_B=\mathbf 0$.
6. The Euler term vanishes for constant angular velocity.
7. A positive maximum heave displacement has downward frame acceleration and upward
   apparent acceleration.
8. Every angle entering `SRFUtils` is in radians.
9. Every cross-product operand is expressed in frame $B$.
10. A source is included exactly once in the momentum equations and in the corresponding
    Rhie-Chow `body_force` configuration.
