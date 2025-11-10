# Conjugate heat transfer over a flat plate

!tag name=Conjugate heat transfer on a flat plate
     image=../../media/validation/free_flow/heat_transfer/flat_plate/icon.png
     description=Conjugate heat transfer on a flat plate
     pairs=flow_type:free-flow
           compressibility:incompressible
           heattransfer:conjugate_heat_transfer
           convection_type:forced
           transient:steady
           flow_regime:laminar
           fluid:air
           number_of_phases:one

## Problem Description

This validation case demonstrates conjugate heat transfer (CHT) between a laminar fluid flow
and a conducting solid wall in a 2D configuration [!cite](luikov1974conjugate). The case is
modeled after the benchmark problem presented in [!cite](verstraete2016stability),
which investigates the stability and accuracy of static coupling approaches for CHT problems.

The problem involves forced convection of air over a solid aluminum wall with a finite
thickness, where heat conduction within the solid is coupled with convective transport in the fluid.
Temperature and heat flux continuity are enforced at the fluid–solid interfaces.

The geometry of the problem is depicted in [!ref](fig:geom).

!style halign=center
!media media/validation/free_flow/heat_transfer/flat_plate/geometry.png style=width:100% id=fig:geom caption=The geometry of the problem (dimensions in m).

We note that an entrance reagion has been added to the domain to allow boundary layer buildup before the
flow reaching the solid domian.
The boundary conditions for the momentum and pressure equations are the following:

!table id=tab:momentum-pressure-bc caption=Momentum and pressure boundary conditions.
| Boundary | Boundary condition |
| --- | --- |
| +Fluid inlet+ | Fixed velocity, $u_{in}=13~\frac{m}{s}$ |
| +Fluid top+ | Free stream |
| +Fluid outlet+      | Outlet pressure, $p_{out}=1.035\times10^5~Pa$ |
| +Fluid bottom and solid top+ | No-slip |

The boundary conditions for the fluid and solid energy equations are the following:

!table id=tab:energy-bc caption=Solid and fluid energy boundary conditions.
| Boundary | Boundary condition |
| --- | --- |
| +Fluid inlet+ | Fixed temperature, $T_{in}=1000~K$ |
| +Fluid top+ | Insulated |
| +Fluid outlet+ | Outlet enthalpy |
| +Fluid bottom before solid+ | Insulated |
| +Solid top+ | Conjugate heat transfer |
| +Solid bottom+ | Fixed temperature $T_b=600~K$ |
| +Solid sides+ | Insulated |

Material properties used are given in [tab:matprops].

!table id=tab:matprops caption=Material properties used in the CHT benchmark case.
| Region | Density, +$\rho$+ $\left(\frac{kg}{m^3}\right)$ | Specific Heat, +$c_p$+ $\left(\frac{J}{kg\cdot K}\right)$ | Conductivity, +$k$+ $\left(\frac{W}{m\cdot K}\right)$ | Dynamic Viscosity, +$\mu$+ $\left(Pa\cdot s\right)$ |
| --- | --- | --- | --- | --- |
| Fluid (air) | 0.3525 | 1142.6 | 0.06808 | 3.95$\times 10^{-5}$ |
| Solid (fictional) | Not used | Not used | 0.2876 | Not used |

## `OpenPronghorn` Model

The mesh is generated using the MOOSE mesh generator system and includes both solid and fluid regions. A
representative section of the mesh is shown in +[fig:mesh]+ below. The cell density near the fluid–solid
interface is refined to capture strong thermal gradients. The refinement is achieved by the introduction of
mesh bias parameters.

!media media/validation/free_flow/heat_transfer/flat_plate/mesh.png style=width:70%;margin-left:auto;margin-right:auto;text-align:center; id=fig:mesh caption=Mesh showing the fluid and solid regions in the 2D CHT channel.

Altogether, the mesh has 125 divisions horizontally in the entrance region and 250 in the fluid region
above the plate. The solid plate has 80 while the fluid domain 160 vertical divisions. Altogether,
the mesh consists of 80,000 quad cells. A mesh sensitivity has been carried out and the included model
provides closely converged results with the least amount of cells.

The simulation employs a +SIMPLE solver+ to obtain the steady state solution. Second order geometric average interpolation has been
used for the advection discretization, while the stress and conduction operators have been
discretized using an arithmetic average discretization. The mesh is orthogonal, so no correction terms were used.
The input file for this case is embedded below.

!listing /validation/free_flow/heat_transfer/flat_plate/cht_rob-rob.i

## Results

Two key validation metrics are examined:

1. **Vertical temperature profile:** temperature variation along the vertical line at $x=0.1~m$.
2. **Wall temperature distribution:** temperature along the solid-fluid interface.

To judge the acceptability of the results, we utilize the analytic expressions by Luikov from
[!cite](luikov1974conjugate). These expressions utilized techniques like the boundary layer integral (BL)
and differential heat transfer (DHT). It is important to note, however, that these expressions were derived
using approximations such as assuming an averaged velocity in the boundary layer and no stream-wise
conduction in the solid and fluid domains.

[!ref](fig:vertprofile) presents the temperature profile along the $x=0.1~m$ vertical line. We see that the
numerical solution is between the two analytic expressions. As expected, the numerical results are slightly
lower than the DHT approach mainly due to the neglect of the stream-wise conduction in the corresponding derivation in DHT.

!media media/validation/free_flow/heat_transfer/flat_plate/vertical_plot.py
       image_name=vertical-temperature.png
       id=fig:vertprofile
       caption=Vertical temperature profile at $x=0.1~m$.
       style=width:70%;margin-left:auto;margin-right:auto;text-align:center

[!ref](fig:interfaceprofile) presents the temperature profile along the solid-fluid interface. We see considerable
differences at the leading edge of the plate where the gradients are the sharpest. Heat conduction in the solid
is especially high serving as an explanation for the experienced deviation.

!media media/validation/free_flow/heat_transfer/flat_plate/interface_plot.py
       image_name=interface-temperature.png
       id=fig:interfaceprofile
       caption=Temperature profile along the fluid-solid interface.
       style=width:70%;margin-left:auto;margin-right:auto;text-align:center

## Discussion of errors

Mainly due to the approximations used in the analytic solutions, we see a considerable
(up to 7-8%) error in the interface temperatures close to the leading edge. Closer to the
middle and end tail of the plate this error shrinks to approximately 2%.

Using the automatic testing of OpenPronghorn, we make sure this error does not increase over time.
We allow a 1% deviation from previous error levels. We emphasize this is a change in the error levels, not in
the solution itself. If software lead to a different solution with errors that surpass this level we
declare a failed validation test. Furthermore, if we encounter changes in results more than 0.1%
we notify the developers with a failed regression test. The validation test might still
pass at this point.

## Performance Chart

!alert note
The following figure showcases the measured total runtime of the problem over the
commit history. We utilized INL's High-Performance Computational resources to run these
simulations so runtimes might vary depending on which physical resource the job got allocated.

!media media/validation/free_flow/heat_transfer/flat_plate/flat_plate_performance.py
       image_name=flat_plate_performance.png
       id=fig:st_performance
       caption=Runtime over the latest commtits.
       style=width:75%;margin-left:auto;margin-right:auto;text-align:center
