# Single Rotating Frame Modeling in OpenPronghorn

This page documents the single rotating frame formulation used in OpenPronghorn for
thermal-hydraulic simulations of systems subjected to prescribed maritime motion.
The formulation can also be extended to rotating systems with a prescribed inclination angle,
angular velocity and acceleration.
The formulation keeps the computational mesh fixed to the system and expresses
the governing equations in that non-inertial body frame. Vessel rotation, translation,
and the changing direction of gravity are represented through momentum sources instead
of mesh motion.

The implementation covers:

- pitch, yaw, and roll transformations between inertial and body coordinates;
- body-frame angular velocity and angular acceleration for coupled rotations;
- centrifugal, Euler, and Coriolis apparent accelerations;
- rotation of gravity and other inertial-frame vector sources;
- Boussinesq buoyancy in the moving frame; and
- harmonic heave as a translational apparent acceleration.

The marine-coordinate and rigid-body notation follows the general conventions described by
[!cite](sname1950nomenclature).

## Reference frames and notation

Two right-handed Cartesian reference frames are used:

- The inertial frame $I$ is fixed in space. Gravity and other environmental vectors are
  normally prescribed in this frame.
- The body frame $B$ is attached to the vessel and reactor. The computational mesh,
  velocity components, and momentum equations are expressed in this frame.

The body axes follow the maritime convention:

| Axis | Direction | Translational motion | Rotational motion |
| :- | :- | :- | :- |
| $x_B$ | longitudinal | surge | roll, $\phi$ |
| $y_B$ | transverse | sway | pitch, $\theta$ |
| $z_B$ | vertical | heave | yaw, $\psi$ |

Let $O_B$ denote the origin of the body frame and motion reference point. In the current
SRF implementation this point is identified with the prescribed metacenter location.
Its inertial position is $\mathbf R_I$. The position of a fluid cell relative to this
point is

\begin{equation}
\mathbf r_B = \mathbf x_B - \mathbf x_{mc,B}.
\label{eq:srf-position-vector}
\end{equation}

The [LinearFVSRFFunctorMaterial.md] object creates this vector as the `r_mc` functor.
All cross products in the SRF acceleration kernels use body-frame components, so
$\mathbf r_B$, $\boldsymbol{\omega}_B$, $\dot{\boldsymbol{\omega}}_B$, and
$\mathbf u_B$ must be expressed in the same frame.

## Pitch-yaw-roll change of axes

### Elementary rotations

The elementary right-handed rotation matrices are

\begin{equation}
\mathbf R_x\phi =
\begin{bmatrix}
1 & 0 & 0 \\
0 & \cos\phi & -\sin\phi \\
0 & \sin\phi & \cos\phi
\end{bmatrix},
\end{equation}

\begin{equation}
\mathbf R_y\theta =
\begin{bmatrix}
\cos\theta & 0 & \sin\theta \\
0 & 1 & 0 \\
-\sin\theta & 0 & \cos\theta
\end{bmatrix},
\end{equation}

and

\begin{equation}
\mathbf R_z\psi =
\begin{bmatrix}
\cos\psi & -\sin\psi & 0 \\
\sin\psi & \cos\psi & 0 \\
0 & 0 & 1
\end{bmatrix}.
\end{equation}

OpenPronghorn uses the body-to-inertial direction-cosine matrix

\begin{equation}
\mathbf C_{IB}
= \mathbf R_y\theta\,\mathbf R_z\psi\,\mathbf R_x\phi.
\label{eq:srf-body-to-inertial}
\end{equation}

Matrix multiplication acts from right to left: a body-frame vector is first acted on by
the roll rotation, then the yaw rotation, and finally the pitch rotation. Rotation order
matters; exchanging any two matrices generally produces a different orientation.

The inertial-to-body transformation is the transpose

\begin{equation}
\mathbf C_{BI}
= \mathbf C_{IB}^{T}
= \mathbf R_x(-\phi)\,\mathbf R_z(-\psi)\,\mathbf R_y(-\theta).
\label{eq:srf-inertial-to-body}
\end{equation}

For an inertial-frame vector $\mathbf v_I$ and the same physical vector expressed in the
body frame,

\begin{equation}
\mathbf v_B = \mathbf C_{BI}\mathbf v_I,
\qquad
\mathbf v_I = \mathbf C_{IB}\mathbf v_B.
\label{eq:srf-vector-transform}
\end{equation}

### Implemented transformation matrix

Using $c_\theta=\cos\theta$, $s_\theta=\sin\theta$,
$c_\psi=\cos\psi$, $s_\psi=\sin\psi$, $c_\phi=\cos\phi$, and
$s_\phi=\sin\phi$, the inertial-to-body matrix implemented in [SRFUtils.md] is

\begin{equation}
\mathbf C_{BI} =
\begin{bmatrix}
c_\psi c_\theta & s_\psi & -c_\psi s_\theta \\
s_\phi s_\theta-c_\phi s_\psi c_\theta & c_\phi c_\psi &
s_\phi c_\theta+c_\phi s_\psi s_\theta \\
c_\phi s_\theta+s_\phi s_\psi c_\theta & -s_\phi c_\psi &
c_\phi c_\theta-s_\phi s_\psi s_\theta
\end{bmatrix}.
\label{eq:srf-implemented-rotation}
\end{equation}

At zero pitch, yaw, and roll, $\mathbf C_{BI}$ reduces to the identity matrix. Because
the matrix is orthogonal, its transpose is also its inverse, and vector magnitudes and
dot products are preserved by the transformation.

## Coupled angular velocity and angular acceleration

Euler-angle rates are not generally the components of angular velocity. For the rotation
sequence in Eq. \eqref{eq:srf-body-to-inertial}, the angular velocity expressed in the
body frame is

\begin{equation}
\boldsymbol{\omega}_B =
\begin{bmatrix}
\dot{\phi} + \sin\psi\,\dot{\theta} \\
\cos\phi\cos\psi\,\dot{\theta} + \sin\phi\,\dot{\psi} \\
-\sin\phi\cos\psi\,\dot{\theta} + \cos\phi\,\dot{\psi}
\end{bmatrix}.
\label{eq:srf-angular-velocity}
\end{equation}

Only when the rotations are uncoupled or the relevant angles vanish can the components
of $\boldsymbol{\omega}_B$ be identified directly with roll, pitch, and yaw rates.

Differentiating Eq. \eqref{eq:srf-angular-velocity} gives the angular acceleration in
the body frame:

\begin{equation}
\dot{\boldsymbol{\omega}}_B =
\begin{bmatrix}
\ddot{\phi}+\sin\psi\,\ddot{\theta}
+\cos\psi\,\dot{\psi}\dot{\theta} \\
\cos\phi\cos\psi\,\ddot{\theta}+\sin\phi\,\ddot{\psi}
-\sin\phi\cos\psi\,\dot{\phi}\dot{\theta}
-\cos\phi\sin\psi\,\dot{\psi}\dot{\theta}
+\cos\phi\,\dot{\phi}\dot{\psi} \\
-\sin\phi\cos\psi\,\ddot{\theta}+\cos\phi\,\ddot{\psi}
-\cos\phi\cos\psi\,\dot{\phi}\dot{\theta}
+\sin\phi\sin\psi\,\dot{\psi}\dot{\theta}
-\sin\phi\,\dot{\phi}\dot{\psi}
\end{bmatrix}.
\label{eq:srf-angular-acceleration}
\end{equation}

These are the expressions used by [LinearFVSRFFunctorMaterial.md] in
`pitch_yaw_roll` mode.

## Harmonic maritime motion

Pitch, yaw, and roll may be prescribed independently as harmonic functions. For
$q\in\{\phi,\theta,\psi\}$,

\begin{equation}
q(t) = A_q\sin\left(\Omega_q t+\delta_q\right),
\qquad
\Omega_q = \frac{2\pi}{T_q},
\label{eq:srf-harmonic-angle}
\end{equation}

\begin{equation}
\dot q(t) = A_q\Omega_q\cos\left(\Omega_q t+\delta_q\right),
\label{eq:srf-harmonic-angle-rate}
\end{equation}

and

\begin{equation}
\ddot q(t) = -A_q\Omega_q^2\sin\left(\Omega_q t+\delta_q\right).
\label{eq:srf-harmonic-angle-acceleration}
\end{equation}

The angle amplitudes and phases accepted by `LinearFVSRFFunctorMaterial` are supplied in
degrees and converted internally to radians. The generated angle functors,
$\boldsymbol{\omega}_B$, and $\dot{\boldsymbol{\omega}}_B$ use radians. Every harmonic
period must be greater than zero.

In `fixed` mode, the orientation angles, angular velocity, and angular acceleration are
specified independently. The user is responsible for ensuring that these quantities are
kinematically consistent when a physically realizable motion history is required.

## Transformation of the momentum equation

### Acceleration transport theorem

The inertial velocity of a material point can be written as

\begin{equation}
\mathbf u_I
= \dot{\mathbf R}_I
+ \mathbf C_{IB}
\left(\mathbf u_B+\boldsymbol{\omega}_B\times\mathbf r_B\right),
\label{eq:srf-velocity-transform}
\end{equation}

where $\mathbf u_B$ is the velocity relative to the body-fixed mesh. Differentiating
Eq. \eqref{eq:srf-velocity-transform} gives the acceleration transport theorem:

\begin{equation}
\mathbf a_I
= \ddot{\mathbf R}_I
+ \mathbf C_{IB}\left[
\frac{D_B\mathbf u_B}{Dt}
+ \dot{\boldsymbol{\omega}}_B\times\mathbf r_B
+ \boldsymbol{\omega}_B\times
  \left(\boldsymbol{\omega}_B\times\mathbf r_B\right)
+ 2\boldsymbol{\omega}_B\times\mathbf u_B
\right].
\label{eq:srf-acceleration-transport}
\end{equation}

The terms inside the brackets represent relative, tangential, centripetal, and Coriolis
acceleration, respectively. Define the translational acceleration of the body-frame
origin expressed in body coordinates as

\begin{equation}
\mathbf a_{O,B}=\mathbf C_{BI}\ddot{\mathbf R}_I.
\label{eq:srf-origin-acceleration}
\end{equation}

Solving Eq. \eqref{eq:srf-acceleration-transport} for the relative acceleration yields

\begin{equation}
\frac{D_B\mathbf u_B}{Dt}
= \mathbf C_{BI}\mathbf a_I
-\mathbf a_{O,B}
-\dot{\boldsymbol{\omega}}_B\times\mathbf r_B
-\boldsymbol{\omega}_B\times
  \left(\boldsymbol{\omega}_B\times\mathbf r_B\right)
-2\boldsymbol{\omega}_B\times\mathbf u_B.
\label{eq:srf-relative-acceleration}
\end{equation}

### Body-frame conservation equations

For a rigid coordinate transformation, volume is preserved and the conservative mass
equation retains its usual form:

\begin{equation}
\frac{\partial\rho}{\partial t}
+\nabla_B\cdot\left(\rho\mathbf u_B\right)=0.
\label{eq:srf-mass}
\end{equation}

The body-frame momentum equation is

\begin{equation}
\frac{\partial\rho\mathbf u_B}{\partial t}
+\nabla_B\cdot\left(\rho\mathbf u_B\otimes\mathbf u_B\right)
=-\nabla_B p
+\nabla_B\cdot\boldsymbol{\tau}_B
+\rho\mathbf g_B
+\mathbf s_B
+\mathbf f_{NI,B},
\label{eq:srf-momentum}
\end{equation}

where

\begin{equation}
\mathbf g_B=\mathbf C_{BI}\mathbf g_I,
\qquad
\mathbf s_B=\mathbf C_{BI}\mathbf s_I,
\label{eq:srf-rotated-physical-sources}
\end{equation}

and the non-inertial force density is

\begin{equation}
\mathbf f_{NI,B}
=-\rho\left[
\mathbf a_{O,B}
+\dot{\boldsymbol{\omega}}_B\times\mathbf r_B
+\boldsymbol{\omega}_B\times
  \left(\boldsymbol{\omega}_B\times\mathbf r_B\right)
+2\boldsymbol{\omega}_B\times\mathbf u_B
\right].
\label{eq:srf-noninertial-force}
\end{equation}

The combined acceleration that acts like an effective gravity field is therefore

\begin{equation}
\mathbf g_{eff,B}
=\mathbf g_B
-\mathbf a_{O,B}
-\dot{\boldsymbol{\omega}}_B\times\mathbf r_B
-\boldsymbol{\omega}_B\times
  \left(\boldsymbol{\omega}_B\times\mathbf r_B\right)
-2\boldsymbol{\omega}_B\times\mathbf u_B.
\label{eq:srf-effective-gravity}
\end{equation}

Because the Coriolis contribution depends on local velocity, $\mathbf g_{eff,B}$ is a
convenient grouping of accelerations rather than a spatially uniform gravitational
field.

## Rotational apparent accelerations

[LinearFVSRFAccelerations.md] implements the rotational part of
Eq. \eqref{eq:srf-noninertial-force}:

\begin{equation}
\mathbf f_{rot,B}
=-\rho\left[
\boldsymbol{\omega}_B\times
  \left(\boldsymbol{\omega}_B\times\mathbf r_B\right)
+\dot{\boldsymbol{\omega}}_B\times\mathbf r_B
+2\boldsymbol{\omega}_B\times\mathbf u_B
\right].
\label{eq:srf-rotational-force}
\end{equation}

The physical meaning and limiting behavior of each term are:

- **Centrifugal apparent force.** The frame-point acceleration
  $\boldsymbol{\omega}_B\times
  \left(\boldsymbol{\omega}_B\times\mathbf r_B\right)$ is centripetal. Moving it to
  the body-frame momentum equation gives the negative, outward centrifugal force.
  It vanishes at the motion origin and scales with distance from that origin.
- **Euler apparent force.** The term
  $-\rho\dot{\boldsymbol{\omega}}_B\times\mathbf r_B$ is caused by changing angular
  velocity. It vanishes for constant $\boldsymbol{\omega}_B$ and at the motion origin.
- **Coriolis apparent force.** The term
  $-2\rho\boldsymbol{\omega}_B\times\mathbf u_B$ acts only on fluid moving relative to
  the body frame. It is perpendicular to both $\boldsymbol{\omega}_B$ and
  $\mathbf u_B$ and therefore does no direct mechanical work in the continuous
  kinetic-energy equation.

The rotational kernel does not include translational acceleration or gravity. Those
effects are introduced by the heave and gravity-related objects described below.

## Gravity in the body frame

Gravity is constant in the inertial frame but generally changes direction in body
coordinates. With

\begin{equation}
\mathbf g_I=
\begin{bmatrix}
0 \\ 0 \\ -g
\end{bmatrix},
\end{equation}

Eq. \eqref{eq:srf-implemented-rotation} gives

\begin{equation}
\mathbf g_B =
\begin{bmatrix}
g\cos\psi\sin\theta \\
-g\left(\sin\phi\cos\theta
+\cos\phi\sin\psi\sin\theta\right) \\
-g\left(\cos\phi\cos\theta
-\sin\phi\sin\psi\sin\theta\right)
\end{bmatrix}.
\label{eq:srf-gravity-components}
\end{equation}

For pure yaw, gravity is unchanged because the rotation is about the inertial vertical
axis. For pure pitch and pure roll, respectively,

\begin{equation}
\mathbf g_B^{pitch}=
\begin{bmatrix}
g\sin\theta \\ 0 \\ -g\cos\theta
\end{bmatrix},
\qquad
\mathbf g_B^{roll}=
\begin{bmatrix}
0 \\ -g\sin\phi \\ -g\cos\phi
\end{bmatrix}.
\label{eq:srf-pure-angle-gravity}
\end{equation}

These limiting forms provide useful checks on sign, angle units, and axis definitions.

## Boussinesq buoyancy under rotation

Under the Boussinesq approximation, density is written as

\begin{equation}
\rho(T)=\rho_{ref}\left[1-\alpha\left(T-T_{ref}\right)\right].
\label{eq:srf-boussinesq-density}
\end{equation}

The gravitational force density becomes

\begin{equation}
\rho\mathbf g_B
=\rho_{ref}\mathbf g_B
-\rho_{ref}\alpha\left(T-T_{ref}\right)\mathbf g_B.
\label{eq:srf-boussinesq-split}
\end{equation}

The first term is the reference hydrostatic contribution and may be absorbed into the
pressure field. [LinearFVSRFMomentumBoussinesq.md] implements the temperature-dependent
part:

\begin{equation}
\mathbf b_{buoy,B}
=-\rho_{ref}\alpha\left(T-T_{ref}\right)\mathbf g_B.
\label{eq:srf-boussinesq-source}
\end{equation}

Gravity is rotated before the selected momentum component is evaluated. A hot fluid with
$T>T_{ref}$ therefore accelerates opposite to $\mathbf g_B$, regardless of the current
vessel orientation.

The heave source is separate from Eq. \eqref{eq:srf-boussinesq-source}. Translational
frame acceleration acts on the complete density and must not be multiplied by
$\alpha\left(T-T_{ref}\right)$.

## Harmonic heave

Heave is a translation of the body-frame origin. For a prescribed displacement along the
selected vertical direction,

\begin{equation}
z_h(t)=A_h\sin\left(\Omega_h t+\delta_h\right),
\qquad
\Omega_h=\frac{2\pi}{T_h},
\label{eq:srf-heave-position}
\end{equation}

the frame acceleration is

\begin{equation}
\ddot z_h(t)
=-A_h\Omega_h^2\sin\left(\Omega_h t+\delta_h\right).
\label{eq:srf-heave-frame-acceleration}
\end{equation}

The apparent acceleration is opposite to the acceleration of the frame:

\begin{equation}
a_h=-\ddot z_h
=A_h\Omega_h^2\sin\left(\Omega_h t+\delta_h\right).
\label{eq:srf-heave-apparent-acceleration}
\end{equation}

[LinearFVSRFMomentumHeave.md] applies

\begin{equation}
b_{h,i}=\rho a_h
\label{eq:srf-heave-source}
\end{equation}

to the momentum component selected by `momentum_component`. With the standard maritime
axes, this is normally the $z_B$ equation. The current heave object applies the scalar
acceleration directly to the selected body component. It does not rotate a translational
acceleration vector through pitch, yaw, and roll.

## Linear finite-volume realization

Each SRF momentum object evaluates a source density in a cell and multiplies it by the
cell volume before adding it to the right-hand side of the corresponding segregated
momentum system. One object is created for each affected momentum component.

| Theory contribution | OpenPronghorn object |
| :- | :- |
| Angles, $\mathbf r_B$, $\boldsymbol{\omega}_B$, and $\dot{\boldsymbol{\omega}}_B$ | [`LinearFVSRFFunctorMaterial.md`](LinearFVSRFFunctorMaterial.md) |
| $\mathbf C_{BI}$ and $\mathbf C_{IB}$ | [SRFUtils.md] |
| Centrifugal, Euler, and Coriolis sources | [`LinearFVSRFAccelerations.md`](LinearFVSRFAccelerations.md) |
| Rotated Boussinesq buoyancy | [`LinearFVSRFMomentumBoussinesq.md`](LinearFVSRFMomentumBoussinesq.md) |
| Harmonic translational heave | [`LinearFVSRFMomentumHeave.md`](LinearFVSRFMomentumHeave.md) |

The source kernels return zero direct matrix contribution. The Coriolis term nevertheless
depends on velocity and is evaluated from the supplied velocity functors, so it is
effectively lagged within the segregated iteration. The SRF kernel names should also be
included in the `body_force` list of the Rhie-Chow interpolator when they contribute to
the pressure-velocity coupling.

## Consistency requirements and common errors

1. **Use one rotation convention.** Equations \eqref{eq:srf-body-to-inertial} through
   \eqref{eq:srf-implemented-rotation} define the implemented order. A standard rotation
   matrix taken from another convention cannot be substituted without also changing the
   angular-velocity transformation.
2. **Use radians internally.** The SRF functor material converts degree-based motion
   inputs to radians. Any independently supplied angle functor must already evaluate in
   radians.
3. **Keep vectors in a common frame.** Every vector in a cross product must use body-frame
   components. Gravity and general vector sources are transformed with $\mathbf C_{BI}$
   before components are extracted.
4. **Use the correct motion origin.** Centrifugal and Euler accelerations depend on
   $\mathbf r_B$. An incorrect `mc_origin` changes both magnitude and direction.
5. **Do not use Euler-angle rates as angular-velocity components.** For coupled motion,
   use Eq. \eqref{eq:srf-angular-velocity}.
6. **Do not combine heave with Boussinesq scaling.** Heave contributes $\rho a_h$ and is
   independent of $T-T_{ref}$ and $\alpha$.
7. **Avoid double counting.** A force represented with `LinearFVSRFSource` should not be
   added again through a dedicated gravity, heave, or acceleration object.
8. **Check limiting cases.** Zero angles must recover $\mathbf C_{BI}=\mathbf I$; zero
   angular velocity removes centrifugal and Coriolis sources; zero angular acceleration
   removes the Euler source; and $\mathbf r_B=\mathbf 0$ removes the centrifugal and
   Euler sources at the motion origin.
