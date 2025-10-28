# Vortex Shedding Over Circle (2D)

!tag name=Vortex Shedding Over Circle
     image=../../media/validation/free_flow/isothermal/vortex_shedding/cylinder/icon.png
     description=Vortex shedding behind a cylinder
     pairs=flow_type:free-flow
                       compressibility:incompressible
                       heattransfer:isothermal
                       convection_type:forced
                       transient:transient
                       flow_regime:laminar
                       fluid:fictional
                       flow_configuration:free-flow
                       number_of_phases:one

## Problem description

This problem describes vortex shedding behind a cylinder in a 2D domain. The
detailed description of the brenchmark can be found in [!cite](schafer1996benchmark).
The problem is slightly asymmetric in the vertical direction to facilitate vortex
shedding. The geometry of the problem is depicted below.

!style halign=center
!media media/validation/free_flow/isothermal/vortex_shedding/cylinder/geometry.png style=width:100% id=fig:geom caption=The geometry of the problem (dimensions in m). The diameter of the cylinder is 0.1 m.

A no-slip boundary condition is used for the top, bottom
and cylinder boundaries. The fluid enters the domain through
the left boundary with a parabolic profile described below:

\begin{equation}
U(y) = \frac{6 (y+0.2)(0.21-y)}{0.41^2}
\end{equation}

An outflow boundary condition is used for the right boundary with fixed dynamic pressure of $0~Pa$.

The transient simulation starts with the fluid at rest.

The fluid itself is fictional with the following material properties:

!table id=tab:matprops caption=Material properties for the benchmark case.
| Parameter | Value |
| --- | --- |
Density | $1~\frac{kg}{m^3}$ |
Dynamic viscosity |  $10^{-3}~Pa\cdot s$ |

## OpenPronghorn Model

The main input parameters for the runs are summarized in a header file
which is presented below:

!listing /validation/free_flow/isothermal/vortex_shedding/cylinder/header.i

The mesh is generated using the native mesh generation capabilities in MOOSE.
The input file is presented below. The mesh contains 21,092 quadriliteral cells altogether and is depicted in Figure REF.

!listing /validation/free_flow/isothermal/vortex_shedding/cylinder/mesh.i

!style halign=center
!media media/validation/free_flow/isothermal/vortex_shedding/cylinder/mesh.png style=width:100% id=fig:mesh caption=The mesh used in this problem.

The input file for the solve is depicted below.

!listing /validation/free_flow/isothermal/vortex_shedding/cylinder/flow.i

As shown, we compute the lift and drag coefficients of the cylinder.
These time series are then analyzed to determine the vortex shedding frequency.

## Results

The main quantity of interest is the Strouhal number, which can be defined as:

\begin{equation}
St = \frac{f D_h}{U}
\end{equation}

where $f$ denotes the frequency, $D_h$ the hydraulic diameter and $U$ the
characteristic speed. The benchmark in [!cite](schafer1996benchmark)
prescribes an acceptable $St$ range of $[0.295, 0.305]$ for this Reynolds number.


!media media/validation/free_flow/isothermal/vortex_shedding/cylinder/strouhal_plot.py
       image_name=strouhal.png
       id=fig:st
       caption=The Strouhal number supplied by OpenPronghorn with respect to the acceptable range.
       style=width:50%;margin-left:auto;margin-right:auto;text-align:center

## Performance Chart

!alert note
The following figure showcases the measured total runtime of the problem over the different
commit history. We utilized INL's High-Performance Computational resources to run these
simulations so runtimes might vary depending on which physical resource the job got allocated.

!media media/validation/free_flow/isothermal/vortex_shedding/cylinder/strouhal_performance.py
       image_name=strouhal_performance.png
       id=fig:st_performance
       caption=Runtime over the latest commits.
       style=width:75%;margin-left:auto;margin-right:auto;text-align:center


