# Turbulent Flow in a 2D Channel

!tag name=Turbulent Flow in a Channel
     image=../../media/validation/free_flow/isothermal/2d_turbulent_channel/1_icon.png
     description=Fully developed, turbulent flow in a 2D channel
     pairs=compressibility:incompressible
           heattransfer:isothermal
           convection_type:forced
           transient:steady
           flow_regime:turbulent
           fluid:fictional
           flow_configuration:free-flow
           number_of_phases:one

## Problem Description

This problem describes a fully-developed, turbulent channel flow in a 2D domain. Two cases were studied based on the following Reynolds numbers, $Re$, and their respective friction Reynolds numbers, $Re_\tau$:

\begin{equation}
Re = 14,000 \; \text{and} \; Re_\tau = 395 \\\\
Re = 22,250 \; \text{and} \; Re_\tau = 590
\end{equation}

## Modeling Parameters

A detailed description of the benchmark can be found in [!cite](kim1987benchmark) or the [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case032).

A no-slip boundary condition is used for the top and bottom wall boundaries. The fluid enters the domain through the left boundary. An outflow boundary condition is used for the right boundary with fixed dynamic pressure of $0~Pa$. The initial horizontal velocity is $1 \times 10^{-8} \frac{m}{s}$ and the initial vertical velocity is $0 \frac{m}{s}$.

The fluid itself is air with the following material properties:

!table id=tab:matprops caption=Material properties for the benchmark case.
| Parameter | Value | Units |
| --- | --- | --- |
Density, +$\rho$+ | $1$ | $\frac{kg}{m^3}$
Dynamic viscosity, +$\mu$+ |  $\frac{\rho vH \cdot 2}{Re}$  |  $Pa \cdot s$  |

+Note:+ The dynamic viscosity varies based on the Reynolds number for each case and was calculated within the input file (see `OpenPronghorn` model).

The k-epsilon turbulence model is implemented using the following initial conditions:

!table id=tab:kepsconds caption= *k-epsilon* turbulence model initial conditions.
| Condition | Equation | Units |
| --- | --- | --- |
Initial turbulence kinetic energy, $k_{init}$ | $\frac{3}{2}(I \cdot |U_{ref}|)^2$ | $\frac{J}{kg}$ |
Initial turbulence dissipation rate, $\varepsilon_{init}$ |  $\frac{C_{\mu}^\frac{3}{4} \cdot k_{init}^\frac{3}{2}}{L}$ | $\frac{J}{kg \cdot s}$ |

where $I$ is the turbulence intensity, $U_{ref}$ is a reference flow speed, and $L$ is a reference length.

## `OpenPronghorn` Model

The mesh is generated using the native mesh generation capabilities in MOOSE. The mesh contains 4,092 quadrilateral cells altogether and is depicted in +[fig:mesh395]+ and +[fig:mesh590]+ below. Notice the mesh for $Re_\tau = 590$ has a different configuration for the intervals in the Y direction. The size of the first cell near the wall is tied to the $y^{+}$ used in the wall function.

!media media/validation/free_flow/isothermal/2d_turbulent_channel/mesh395_2.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh395 caption=Closeup of mesh at the outlet for $Re_\tau = 395$ (channel half-width of 1~m).

!media media/validation/free_flow/isothermal/2d_turbulent_channel/mesh590_3.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh590 caption=Closeup of mesh at the outlet for $Re_\tau = 590$ (channel half-width of 1~m).

The simulations were executed using a +linear SIMPLE finite volume solver+ with +k-epsilon turbulence modeling+.
Additional model and discretization parameters may be found in [tab:numparameters].

!table id=tab:numparameters caption=Model and discretization parameters.
| Parameter | Inlet | Outlet | Walls | Bulk Face Interpolation |
| --- | --- | --- | --- | --- |
| +Velocity+ | Dirichlet          | Fully-developed | No-slip | Rhie-Chow mass flux, upwind advected quantity |
| +Pressure+ | Two-term expansion | Dirichlet | One-term expansion | Average |
| +TKE+      | Dirichlet          | Fully-developed | Equilibrium wall function | Upwind |
| +TKED+.    | Dirichlet          | Fully-developed | Equilibrium wall function | Upwind |
| $\mu_t$     |                    |                 | Equilibrium wall function | Average  |

The input file for the solve with $Re_\tau=395$ is embedded below.

!listing validation/free_flow/isothermal/ercoftac_032_channel_flow/Ret395_linear_SIMPLE_k-epsilon.i

## Results

The main quantities of interest are the horizontal velocities, $U_{x}$, at the channel outlet and the $y^{+}$.
The $y^{+}$ in the first cell is used to compute the normalization factor for the velocity.
The horizontal velocities were recorded at varying heights up the channel and normalized by the frictional velocity $u_\tau$ via the following operations:

\begin{equation}
u_\tau = \frac{y^{+} \cdot Re_{bulk}}{\Delta y}
\end{equation}

\begin{equation}
U_{norm} = \frac{U_x}{u_\tau}
\end{equation}

where $Re_{bulk}$ is the bulk Reynolds number, and $\Delta y$ is the centroid distance from the wall. The normalized velocities were plotted against the varying heights of the channel and the $y^{+}$ values to compare the `OpenPronghorn` simulation results against the benchmark results.

!media media/validation/free_flow/isothermal/2d_turbulent_channel/channel_plot.py
       image_name=7_Ret395_dual_plots.png
       id=fig:mesh395-2 style=width:80%;margin-left:auto;margin-right:auto;text-align:center
       caption=Axial velocity radial profiles (left) and Law of the wall (right) for $Re_t = 395$.

!media media/validation/free_flow/isothermal/2d_turbulent_channel/channel_plot.py
       image_name=8_Ret590_dual_plots.png id=fig:mesh590-2
       style=width:80%;margin-left:auto;margin-right:auto;text-align:center
       caption=Axial velocity radial profiles (left) and Law of the wall (right) for $Re_t = 590$.

## Validation

The `OpenPronghorn` results are examined using the difference to the normalized axial velocity, calculated as follows:

\begin{equation}
U_{DNS_i}-U_{norm_i}
\end{equation}

where +$U_{DNS_i}$+ is the tangeantial velocity of the benchmark case, +$U_{norm_i}$+ is the velocity of the `OpenPronghorn` model.
The allowable interval is within 2\% of the current solution difference.

The $y^{+}$ of the first cell near the wall at the outlet of the system, in the fully developed flow region is also tracked.
A 2\% deviation is allowed on that value.
