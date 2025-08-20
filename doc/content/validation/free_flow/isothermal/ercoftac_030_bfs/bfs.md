# Backward-Facing Step with Inclined Opposite Wall (2D)

!tag name=BFS with Inclined Opposite Wall
    image=media/validation/free_flow/isothermal/title_image.png
    description=Backward-Facing Wall with Inclined Opposite Wall
    pairs=flow_type:free-flow
                       compressibility:incompressible
                       heattransfer:isothermal
                       convection_type:forced
                       transient:transient
                       flow_regime:turbulent
                       fluid:air
                       flow_configuration:free-flow
                       number_of_phases:one

## Problem Description

This problem describes a fully-developed turbulent flow in a channel with a rear-facing step in a 2D domain. A detailed description of the benchmark can be found in [!cite](driver1985benchmark) or the [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case030). The original benchmark by Driver and Seegmiller tested the case at varying degrees of inclination of the opposite wall. The `OpenPronghorn` model focuses only on the scenario in which the opposing wall is at $0^{\circ}$ (completely horizontal).

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/title_image.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh3 caption=TKED profile near the step.

## Modeling Parameters

A no-slip boundary condition is used for all walls. The fluid enters the domain through the left boundary (inlet) and exits out the right boundary (outlet). The initial conditions including horizontal and vertical velocities, pressure, $k_{init}$, and $\varepsilon_{init}$ are initialized from `.csv` files.

The fluid is air with the following material properties:

!table id=tab:matprops caption=Material properties for the benchmark case.
| Parameter                | Value                   | Units            |
| ---                      | ---                     | ---              |
| Density, $\rho$          | $1.18415$               | $\frac{kg}{m^3}$ |
| Dynamic viscosity, $\mu$ | $1.8551 \times 10^{-5}$ | $Pa \cdot s$     |

The k-epsilon turbulence model is implemented using the following parameters and initial conditions:

!table id=tab:kepsparam caption=*k-epsilon* turbulence model parameters.
| Parameter                         | Value    | Units         |
| ---                               | ---      | ---           |
| Turbulence intensity, $I$         | $0.01$   | Dimensionless |
| Bulk velocity at inlet, $U_{ref}$ | $44.2$   | $\frac{m}{s}$ |
| Reference length, $D$             | $0.2032$ | $m$           |

!table id=tab:kepsconds caption= *k-epsilon* turbulence model initial conditions.
| Condition                                               | Equation                                                   | Units                  |
| ---                                                     | ---                                                        | ---                    |
Initial turbulence kinetic energy, $k_{init}$             | $\frac{3}{2}(I \cdot |U_{ref}|)^2$                         | $\frac{J}{kg}$         |
Initial turbulence dissipation rate, $\varepsilon_{init}$ | $\frac{C_{\mu}^\frac{3}{4} \cdot k_{init}^\frac{3}{2}}{D}$ | $\frac{J}{kg \cdot s}$ |

## `OpenPronghorn` Model

The mesh is loaded externally via `FileMeshGenerator` and contains 20,022 quadrilateral cells altogether. [fig:mesh1] shows the mesh in its entirety. [fig:mesh2] is a closeup of the step and the surrounding mesh profile.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/mesh.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh1 caption=Overview of mesh.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_closeup.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh2 caption=Closeup of the step.

The simulations were executed using a +nonlinear SIMPLE finite volume solver+ with +k-epsilon turbulence modeling+. Additional modeling and discretization parameters are found in [tab:numparameters].

!table id=tab:numparameters caption=Modeling and discretization parameters.
| Parameter  | Inlet              | Outlet          | Walls           | Bulk Face Interpolation             |
| ---        | ---                | ---             | ---             | ---                                 |
| +Velocity+ | Dirichlet          | Fully-developed | No-slip         | Rhie-Chow velocity, Upwind momentum | 
| +Pressure+ | Two-term expansion | Dirichlet       | -               | Average                             |
| +TKE+      | Dirichlet          | Fully-developed | Non-equilibrium | Upwind                              |
| +TKED+.    | Dirichlet          | Fully-developed | Non-equilibrium | Upwind                              |

The input file for the solve is embedded below.

!listing validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input.i

## Results

The main quantities of interest are pressure coefficient $c_{p}$, wall-skin friction coefficient $c_{f}$, and the x-velocity profiles up the height of the channel. 
Further manipulation in post was required to derive $c_{p}$ and $c_{f}$ as `OpenPronghorn`'s ability to directly calculate $c_{p}$ and $c_{f}$ is still in development.

To derive the pressure coefficient $c_{p}$, the pressure variable was recorded by `OpenPronghorn` and divided by the dynamic pressure $q$, as follows:

\begin{equation}
q = \frac{\rho u^2}{2}
\end{equation}

\begin{equation}
c_p = \frac{p_{s}}{q}
\end{equation}

!table id=tab:cpparam caption= $c_{p}$ calculation variables.
| Variable           | Value     | Units            |
| ---                | ---       | ---              |
| Density, $\rho$    | $1.18415$ | $\frac{kg}{m^3}$ |
| Flow velocity, $u$ | $48.18$   | $\frac{m}{s}$    |

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/plots_cp_2.png style=width:50%;margin-left:auto;margin-right:auto;text-align:center; id=fig:plot1 caption=Graph of the pressure coefficient across the length of the channel.

To derive the wall-skin friction coefficient $c_{f}$, the following variables were recorded in `OpenPronghorn` and manipulated as follows: turbulent dynamic viscosity $\mu_{t}$, wall distance $d_{wall}$, and horizontal velocity $u_x$. 

\begin{equation}
q = \frac{\rho u^2}{2}
\end{equation}

\begin{equation}
c_f = \frac{\mu_{t}u_x}{d_{wall}} \cdot \frac{1}{q}
\end{equation}

!table id=tab:cfparam caption= $c_{f}$ calculation variables.
| Variable           | Value     | Units            |
| ---                | ---       | ---              |
| Density, $\rho$    | $1.18415$ | $\frac{kg}{m^3}$ |
| Flow velocity, $u$ | $47$      | $\frac{m}{s}$    |

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/plots_cf_2.png style=width:50%;margin-left:auto;margin-right:auto;text-align:center; id=fig:plot2 caption=Graph of the wall-skin friction coefficient across the length of the channel.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/plots_u_profiles_2.png style=width:100%;margin-left:auto;margin-right:auto;text-align:center; id=fig:plot3 caption=Graphs of the vertical x-velocity profiles at different lengths down the channel.

## Validation

The `OpenPronghorn` results were validated by comparing against the benchmark data using Root Mean Square Error (RMSE) and error bars. The RMSE is calculated as follows:

\begin{equation}
RMSE = \sqrt\frac{\sum(X_{exp_i}-X_{sim_i})^2}{N}
\end{equation}

where +$X_{exp_i}$+ is the benchmark data, +$X_{sim_i}$+ is the `OpenPronghorn` model data, and +$N$+ is the number of data points.

!table id=tab:rmse caption=RMSE for `OpenPronghorn` simulations.
| Parameter                      | RMSE  |
| ---                            | ---   |
| Pressure Coefficient           | 0.220 |
| Wall-Skin Friction Coefficient | 0.145 |
| Velocity at x/H = 1            | 
| Velocity at x/H = 4            | 
| Velocity at x/H = 6            | 
| Velocity at x/H = 10           | 