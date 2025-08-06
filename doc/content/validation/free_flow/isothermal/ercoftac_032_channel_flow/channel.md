# Turbulent Flow in a 2D Channel

!tag name=Turbulent Flow in a Channel
     image=../../media/validation/free_flow/isothermal/2d_turbulent_channel/1_icon.png
     description=Fully developed, turbulent flow in a 2D channel 
     pairs=flow_type:free-flow
                       compressibility:incompressible
                       heattransfer:isothermal
                       convection_type:forced
                       transient:transient
                       flow_regime:turbulent
                       fluid:fictional
                       flow_configuration:free-flow
                       number_of_phases:one

## Problem Description

This problem describes a fully-developed, turbulent channel flow in a 2D domain. Two cases were studied based on the following Reynolds numbers, $Re$, and their respective friction Reynolds numbers, $Re_\tau$: 

\begin{equation} 
Re = 14,000 \; and \; Re_\tau = 395 \\\\
Re = 22,250 \; and \; Re_\tau = 590
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

+Note:+ The dynamic viscosity varies based on the Reynolds number for each case and was calculated within the input file (see `OpenPronghorn` Model).

The k-epsilon turbulence model is implemented using the following initial conditions:

!table id=tab:kepsconds caption= *k-epsilon* turbulence model initial conditions.
| Condition | Equation | Units |
| --- | --- | --- |
Initial turbulence kinetic energy, $k_{init}$ | $\frac{3}{2}(I \cdot |U_{ref}|)^2$ | $\frac{J}{kg}$ |
Initial turbulence dissipation rate, $\varepsilon_{init}$ |  $\frac{C_{\mu}^\frac{3}{4} \cdot k_{init}^\frac{3}{2}}{L}$ | $\frac{J}{kg \cdot s}$ |

where $I$ is the turbulence intensity, $U_{ref}$ is a reference flow speed, and $L$ is a reference length.

## `OpenPronghorn` Model

The mesh is generated using the native mesh generation capabilities in MOOSE. The mesh contains 4,092 quadrilateral cells altogether and is depicted in +Figure 1+ and +Figure 2+ below. Notice the mesh for $Re_\tau = 590$ has a different configuration for the intervals in the Y direction.

!media media/validation/free_flow/isothermal/2d_turbulent_channel/mesh395_2.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh395 caption=Closeup of mesh at the outlet for $Re_\tau = 395$ (dimensions given as wall units). 

!media media/validation/free_flow/isothermal/2d_turbulent_channel/mesh590_3.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh590 caption=Closeup of mesh at the outlet for $Re_\tau = 590$ (dimensions given as wall units).

The simulations were executed using a +linear SIMPLE finite volume solver+ with +k-epsilon turbulence modeling+.
Additional model and discretization parameters may be found in [tab:numparameters].

!table id=tab:numparameters caption=Model and discretization parameters.
| Parameter | Inlet | Outlet | Walls | Bulk Face Interpolation |
| --- | --- | --- | --- | --- |
| +Velocity+ | Dirichlet          | Fully-developed | No-slip | Rhie-Chow velocity, upwind momentum | 
| +Pressure+ | Two-term expansion | Dirichlet | - | Average |
| +TKE+      | Dirichlet          | Fully-developed | Equilibrium | upwind |
| +TKED+.    | Dirichlet          | Fully-developed | Equilibrium | upwind |

The input file for the solve is embedded below.

!listing validation/free_flow/isothermal/ercoftac_032_channel_flow/1_input/Ret395_linear_SIMPLE_k-epsilon.i

## Results

The main quantities of interest are the horizontal velocities, $U_{x}$, at the channel outlet and the $y^{+}$. The horizontal velocities were recorded at varying heights up the channel and normalized by the frictional velocity $u_\tau$ via the following operations: 

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

The `OpenPronghorn` results were validated using the root mean square error, +$RMSE$+, calculated as follows:

\begin{equation}
RMSE = \sqrt\frac{\sum(U_{DNS_i}-U_{norm_i})^2}{N}
\end{equation}

where +$U_{DNS_i}$+ is the velocity of the benchmark case, +$U_{norm_i}$+ is the velocity of the `OpenPronghorn` model, and +$N$+ is the number of data points. The lower and upper validation bounds for the +$RMSE$+ are 0.0 and 0.25, respectively, and assigned based on an average deviation within 2.5% of the ERCOFTAC results.

The RMSE results for the two cases are as follows:

!table id=tab:rmse caption=RMSE for `OpenPronghorn` simulations.
| Re@t@ | RMSE | Acceptable Range |
| --- | --- | --- |
395 | 0.220 | 0.00 -- 0.250 |
590 | 0.145 | 0.00 -- 0.250 |