# Aerosol Particle Distribution and Deposition in a Ventilated Chamber

!tag name=Ventilated Chamber Aerosol Deposition
     image=../../media/validation/aerosols/icon.png
     description=Aerosol transport with gravitational settling and near-wall deposition in a ventilated chamber
     pairs=compressibility:incompressible
           heattransfer:isothermal
           convection_type:forced
           transient:steady
           flow_regime:turbulent
           fluid:air
           flow_configuration:free-flow
           number_of_phases:two

## Scope and Purpose

This validation test exercises `OpenPronghorn`’s capability to predict:

- A steady, turbulent ventilation flow field in an enclosed chamber using an RNG RANS *k-ε* closure with wall functions.

- Eulerian aerosol transport using a drift–flux formulation that includes:

  - advection by the carrier flow,
  - gravitational settling (drift velocity),
  - Brownian + turbulent diffusion, and
  - near-wall deposition represented as a boundary flux condition.

The test is intended to validate the coupled behavior that matters for indoor/aerosol dispersion problems:

1. the resolved mean flow field that transports particles through the enclosure and
2. the near-wall mass transfer / deposition behavior that controls concentration gradients and wall losses.

*Example application*: transport of radioisotopes in aerosol form in the containment of a nuclear reactor.

*Out of scope in this test*: thermal buoyancy, humidity effects, coagulation, resuspension, and non-spherical particle effects.

## Problem Description

Chen et al. studied particle distribution and deposition in a mechanically ventilated “scale room” [!cite](chen2006driftflux).
Air flow loaded with aerosol particles enters via a top window in the left of the domain as shown in +[fig:vel_field]+.
Then, the flow recirculates in the cavity and exits the domain via a window places in the bottom of the right wall in +[fig:vel_field]+.
Aerosol recirculates in the cavity and deposits mainly in the floor of the cavity via gravitational decantation as shown in +[fig:aerosol_field]+.
Also, some of the aerosol deposits at the walls and roof of the domain.

!media media/validation/aerosols/velocity_field.png
       style=width:50%;margin-left:auto;margin-right:auto;text-align:center; id=fig:vel_field
       caption=Velocity field in ventilated cavity.

!media media/validation/aerosols/aerosol_field.png
       style=width:50%;margin-left:auto;margin-right:auto;text-align:center; id=fig:aerosol_field
       caption=Iso-contour plot of aerosol concentration field in ventilated cavity.

### Geometry and Operating Conditions (benchmark)

The chamber geometry is a rectangular enclosure with a square inlet and outlet:

!table id=tab:chen-geom caption=Benchmark geometry and measurement locations. [!cite](chen2006driftflux)
| Item                | Value                                               |
| ------------------- | --------------------------------------------------- |
| Chamber dimensions  | $L \times W \times H = 0.8 \times 0.4 \times 0.4$ m |
| Inlet / outlet size | $0.04 \times 0.04$ m                                |
| Inlet center        | $(x,y,z)=(0.0,\ 0.2,\ 0.36)$ m                      |
| Outlet center       | $(x,y,z)=(0.8,\ 0.2,\ 0.04)$ m                      |
| Center plane        | $y = 0.2$ m                                         |
| Profile locations   | $x = 0.2,\ 0.4,\ 0.6$ m on center plane             |

Chen et al. examined multiple inlet flow rates; this `OpenPronghorn` validation focuses on the **lower inlet velocity** case,
$U_{in}=0.225~\mathrm{m/s}$, which is reported by Chen et al. as corresponding to 10 air changes per hour [!cite](chen2006driftflux).

## `OpenPronghorn` Model

This test uses a steady, incompressible finite-volume formulation with a SIMPLE-based segregated pressure–velocity solve and an RNG RANS *k-ε* turbulence model.
Aerosols are modeled as an Eulerian scalar field with drift–flux transport and a wall deposition boundary condition.
The user can refer to the theory section of aerosol modeling for more details.

### Governing equations and closures

The `OpenPronghorn` model solves:

- Incompressible RANS momentum + continuity for the carrier phase (air),
- Transport equations for turbulent kinetic energy $k$ and dissipation rate $\varepsilon$,
- Aerosol drift–flux advection–diffusion with gravitational settling and wall deposition.

The aerosol effective diffusivity is modeled as:

\begin{equation}
D_\mathrm{eff} = D_B + \frac{\mu_t}{\rho\, Sc_{t,p}},
\end{equation}

where $D_B$ is Brownian diffusivity, $\mu_t$ is turbulent viscosity from the RANS closure, $\rho$ is air density, and $Sc_{t,p}$ is the turbulent Schmidt number for particles (set to 1.0 in this case).

The input file for this validation test is embeded below.

!listing validation/aerosols/chen_steady.i

### Boundary and initial conditions

- **Inlet:** fixed velocity ($U_{in}$) at the inlet patch; turbulence quantities are prescribed from an inlet turbulence intensity and length scale; aerosol concentration is set to $C/C_{in}=1$.
- **Outlet:** outflow condition for velocity and aerosol.
- **Walls:** no-slip with wall functions for turbulence; aerosol wall deposition is imposed via `LinearFVAerosolDepositionBC`.

### Modeled particle properties

The validation targets the same particle size used for the profile comparisons in Chen et al. for the 10$\mu m$ class [!cite](chen2006driftflux).

!table id=tab:chen-particles caption=Particle properties used in the `OpenPronghorn` input.
| Property                 |     Symbol |                             Value |
| ------------------------ | ---------: | --------------------------------: |
| Particle diameter        |      $d_p$ | $1.0\times10^{-5}$ m (10 $\mu m$) |
| Particle density         |   $\rho_p$ |                     1400 kg/m$^3$ |
| Mean free path (air)     |  $\lambda$ |              $6.6\times10^{-8}$ m |
| Turbulent Schmidt number | $Sc_{t,p}$ |                               1.0 |

## Results

The quantities of interest are:

1. *Axial velocity profiles* $u/U_{in}$ vs. height on the chamber center plane at $x=0.2, 0.4, 0.6$ m.
2. *Aerosol concentration profiles* $C/C_{in}$ vs. height at the same locations.

The plots below compare:

- digitized experimental profiles, and
- the corresponding `OpenPronghorn` sampler profiles from the current run.

!media media/validation/aerosols/chen_plot.py
       image_name=chen_profiles.png
       style=width:95%;margin-left:auto;margin-right:auto;text-align:center
       id=fig:chen
       caption=Axial velocity (top row) and aerosol concentration (bottom row) profiles on the center plane.

The mesh used in the validation test is fine enough to produce reasonably accurate results,
yet coarse enough to allow the test to run within a reasonable time frame.
The axial velocity profile is captured with good accuracy, although slight over-diffusion occurs downstream due to the coarse mesh.
The aerosol concentration is accurate on average, indicating that aerosol outflow, along with
deposition on the floors, walls, and roof, is correctly represented.
Furthermore, the concentration profiles align with the experimental results within the bounds of experimental uncertainty.

## Validation

This test follows the same “error-to-experiment envelope” pattern used by the ERCOFTAC validations:

1. Interpolate digitized experimental data onto the sampler heights.
2. Compute the **current error-to-experiment** at each height for velocity and concentration:
   \begin{equation}
   e_\mathrm{cur}(z/H) = \left|q_\mathrm{cur}(z/H) - q_\mathrm{exp}(z/H)\right|
   \end{equation}
3. Compute the **reference (gold) error-to-experiment**:
   \begin{equation}
   e_\mathrm{ref}(z/H) = \left|q_\mathrm{ref}(z/H) - q_\mathrm{exp}(z/H)\right|
   \end{equation}
4. Require the current error to remain within a tolerance band around the reference error:
   \begin{equation}
   (1+\ell)\, e_\mathrm{ref} \le e_\mathrm{cur} \le (1+u)\, e_\mathrm{ref}
   \end{equation}

For this validation, the default bounds are:

- $\ell = -0.02$
- $u = +0.02$

i.e., the current error-to-experiment must stay within **±2%** of the reference error-to-experiment at each sampled height, for each profile location and each quantity (velocity and concentration).

The checks are implemented in `chen_profiles_vnv.py` and configured in the local `validation` file within the test directory.
