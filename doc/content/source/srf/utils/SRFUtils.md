# SRFUtils

## Description

`SRFUtils` provides vector transformations between the inertial metacenter frame and the
body reference frame. The utilities are used by the SRF momentum kernels whenever a
vector defined in one frame must be expressed in the other.

## Rotation convention

Let $\theta$, $\psi$, and $\phi$ denote pitch, yaw, and roll, respectively, and define
$c_p=\cos(\theta)$, $s_p=\sin(\theta)$, $c_y=\cos(\psi)$,
$s_y=\sin(\psi)$, $c_r=\cos(\phi)$, and $s_r=\sin(\phi)$. The inertial-to-body
transformation implemented by `rotateVectorInertialToBody` is

\begin{equation}
\boldsymbol{R}_{I\rightarrow B} =
\begin{bmatrix}
c_yc_p & s_y & -c_ys_p \\
s_rs_p-c_rs_yc_p & c_rc_y & s_rc_p+c_rs_ys_p \\
c_rs_p+s_rs_yc_p & -s_rc_y & c_rc_p-s_rs_ys_p
\end{bmatrix}.
\end{equation}

For a vector $\boldsymbol{v}_I$ expressed in the inertial frame,

\begin{equation}
\boldsymbol{v}_B = \boldsymbol{R}_{I\rightarrow B}\boldsymbol{v}_I.
\end{equation}

The inverse transformation is implemented by `rotateVectorBodyToInertial` and uses the
transpose of the same orthogonal matrix:

\begin{equation}
\boldsymbol{v}_I = \boldsymbol{R}_{I\rightarrow B}^{T}\boldsymbol{v}_B.
\end{equation}

All three angles passed to these functions must be in radians.

## Available functions

```cpp
RealVectorValue rotateVectorInertialToBody(const RealVectorValue & vector,
                                           const Real & pitch_angle,
                                           const Real & yaw_angle,
                                           const Real & roll_angle);

RealVectorValue rotateVectorBodyToInertial(const RealVectorValue & vector,
                                           const Real & pitch_angle,
                                           const Real & yaw_angle,
                                           const Real & roll_angle);
```

Both functions are defined in the `NS::SRF` namespace.

## Example

```cpp
const RealVectorValue gravity_brf =
    NS::SRF::rotateVectorInertialToBody(gravity, pitch_angle, yaw_angle, roll_angle);
```
