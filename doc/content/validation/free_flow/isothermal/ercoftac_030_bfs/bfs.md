# Backward-Facing Step with Straight Opposite Wall (2D)

!tag name=BFS with Inclined Opposite Wall
    image=../../media/validation/free_flow/isothermal/ercoftac_030_bfs/vel_streamlines.png
    description=Backward-Facing Wall with Inclined Opposite Wall
    pairs=flow_type:free-flow
          compressibility:incompressible
          heattransfer:isothermal
          convection_type:forced
          transient:steady
          flow_regime:turbulent
          fluid:air
          flow_configuration:free-flow
          number_of_phases:one


## Problem Description

This problem describes a fully-developed turbulent flow in a channel with a rear-facing step in a 2D domain. A detailed description of the benchmark can be found in [!cite](driver1985benchmark) or the [ERCOFTAC database](http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case030). The original benchmark by Driver and Seegmiller tested the case at varying degrees of inclination of the opposite wall. The `OpenPronghorn` model focuses only on the scenario in which the opposing wall is at $0^{\circ}$ (completely horizontal).

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/geometry.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:geom caption=Backwards-facing step geometry.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/vel_streamlines.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh3 caption=Velocity streamlines near the step.

## Modeling Parameters

The fluid is air with the following material properties:

!table id=tab:matprops caption=Material properties for the benchmark case.
| Parameter                | Value                   | Units            |
| ---                      | ---                     | ---              |
| Density, $\rho$          | $1.18415$               | $\frac{kg}{m^3}$ |
| Dynamic viscosity, $\mu$ | $1.8551 \times 10^{-5}$ | $Pa \cdot s$     |

The standard k-epsilon turbulence model is used with the following parameters and initial conditions:

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


A no-slip boundary condition is used for all walls. The fluid enters the domain through the left boundary (inlet) and exits out the right boundary (outlet).
The inlet velocity and turbulence field profiles have been previously obtained from a well-developed channel flow simulation and are imposed at the inlet boundary using `CSV` files.

!alert note title=Regenerating inlet conditions
An additional input file for a flow channel with the mesh discretization as the inlet of the backward-facing step simulation can be found in the input file.
The `CSV` output of the outlet sampler vectorpostprocessor can be used to generate a new `CSV` input for boundary conditions. Performing this step once
with the flow channel simulation (=recycling the boundary condition) helps develop the flow more before using it at the backward-facing step inlet.

## `OpenPronghorn` Model

The mesh is loaded externally via `FileMeshGenerator` and contains 20,022 quadrilateral cells altogether. +[fig:mesh1]+ shows the mesh in its entirety. +[fig:mesh2]+ is a closeup of the step and the surrounding mesh profile.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/mesh.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh1 caption=Overview of the entire mesh.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_closeup.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh2 caption=Closeup of the step.

The simulations were executed using the +SIMPLE finite volume solver+ with the +standard two-equation k-epsilon turbulence model+. Additional modeling and discretization parameters are found in [tab:numparameters] .

!table id=tab:numparameters caption=Modeling and discretization parameters.
| Parameter  | Inlet              | Outlet          | Walls                         | Bulk Face Interpolation             |
| ---        | ---                | ---             | ---                           | ---                                 |
| +Velocity+ | Dirichlet          | Fully-developed | No-slip                       | Rhie-Chow mass flux, Upwind advected quantity |
| +Pressure+ | Two-term expansion | Dirichlet       | One-term expansion            | Average                             |
| +TKE+      | Dirichlet          | Fully-developed | Non-equilibrium wall function | Upwind                              |
| +TKED+.    | Dirichlet          | Fully-developed | Non-equilibrium wall function | Upwind                              |
| $\mu_t$    |                    |                 | Non-equilibrium wall function | Average                             |


The input file for the solve is embedded below.

!listing validation/free_flow/isothermal/ercoftac_030_bfs/bfs_input.i

## Results

The main quantities of interest are pressure coefficient $c_{p}$, wall skin friction coefficient $c_{f}$, and the x-velocity profiles up the height of the channel.
Further computations using the postprocessor values was required to derive $c_{p}$ and $c_{f}$.

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
| Flow velocity, $u$ | $44.02$   | $\frac{m}{s}$    |

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_plot.py
       image_name=plots_cp_main.png
       style=width:50%;margin-left:auto;margin-right:auto;text-align:center;
       id=fig:plot1
       caption=Graph of the pressure coefficient across the length of the channel.

To derive the wall skin friction coefficient $c_{f}$, the following variables were recorded in `OpenPronghorn` and manipulated as follows: turbulent dynamic viscosity $\mu_{t}$, $\mu_{eff} = \mu + \mu_{t}$ the effective dynamic viscosity, wall distance $d_{wall}$, and horizontal velocity $u_x$. +[fig:plot3]+ shows the horizontal velocity at different lengths along the channel,
and +[fig:plot2]+ shows the wall skin friction coefficient.

\begin{equation}
q = \frac{\rho u^2}{2}
\end{equation}

\begin{equation}
c_f = \frac{\mu_{wall} u_x}{d_{wall}} \cdot \frac{1}{q}
\end{equation}

!table id=tab:cfparam caption= $c_{f}$ calculation variables.
| Variable           | Value     | Units            |
| ---                | ---       | ---              |
| Density, $\rho$    | $1.18415$ | $\frac{kg}{m^3}$ |
| Flow velocity, $u$ | $44.02$      | $\frac{m}{s}$    |

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_plot.py
       image_name=plots_cf_main.png
       style=width:50%;margin-left:auto;margin-right:auto;text-align:center
       id=fig:plot2
       caption=Graph of the wall skin friction coefficient across the length of the channel.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_plot.py
       image_name=plots_u_profiles_main.png
       style=width:100%;margin-left:auto;margin-right:auto;text-align:center
       id=fig:plot3
       caption=Graphs of the horizontal x-velocity profiles at different lengths down the channel.

!alert note
The secondary recirculation zone at the corner of the step is not captured by the standard k-epsilon model with high y$^+$ wall functions.

## Validation

Current `OpenPronghorn` results were validated by comparing against errors between the benchmark data and reference `OpenPronghorn` results. The measurement uncertainty was also
considered, and compared to the error. The measurement uncertainty can be found in [fig:plot2], which can be compared to the errors found in [fig:plot4].
The limited accuracy obtained by the $k-\epsilon$ model was confirmed using external software. The $k-\omega$ model, not yet released in MOOSE,
has been found to have more accurate results.

The errors on this validation case should not increase by more than 1% over the current difference to the validation data.
If the current `OpenPronghorn` simulation data fall within the error ranges, the results are validated. This
ensures the error is only improved over time with refinements in the model. The logic for ensuring the error does not increase
is reproduced below. We first compute the maximum absolute deviation with the 1% allowed increase, then compare
the current error to the maximum absolute deviation. This ensures the performance of the validation basis is not degraded
when modifying `OpenPronghorn` or `MOOSE`.

\begin{equation}
\text{Maximum Absolute Deviation} = 1.01 \cdot |X_R - X_E|
\end{equation}

\begin{equation}
\text{Acceptable Range} = X_E \: \pm \: \text{Maximum Absolute Deviation}
\end{equation}

where $X_E$ is the ERCOFTAC data and $X_R$ is the reference data.
+[fig:plot4]+, +[fig:plot5]+ and +[fig:plot6]+ show the reference errors for the metrics tracked. The current errors are checked
by the validation suite, but are not plotted on this page.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_plot.py
       image_name=plots_cp_cf_error_main.png
       style=width:100%;margin-left:auto;margin-right:auto;text-align:center
       id=fig:plot4
       caption=Reference `OpenPronghorn` error for $c_p$ and $c_f$.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_plot.py
       image_name=plots_u_profiles_error1_main.png
       style=width:100%;margin-left:auto;margin-right:auto;text-align:center
       id=fig:plot5
       caption=Reference `OpenPronghorn` error for vertical x-velocities at x/H = 1 and 4.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_plot.py
       image_name=plots_u_profiles_error2_main.png
       style=width:100%;margin-left:auto;margin-right:auto;text-align:center
       id=fig:plot6
       caption=Reference `OpenPronghorn` error for vertical x-velocities at x/H = 6 and 10.

## Performance Chart

!alert note
The following figure showcases the measured total runtime of the problem over the different
commit history. We utilized INL's High-Performance Computational resources to run these
simulations so runtimes might vary depending on which physical resource the job got allocated.

!media media/validation/free_flow/isothermal/ercoftac_030_bfs/bfs_performance.py
       image_name=bfs_performance.png
       id=fig:bfs_performance
       caption=Runtime over the latest commits.
       style=width:75%;margin-left:auto;margin-right:auto;text-align:center

