# Particle Distribution and Deposition in a Ventilated Chamber (Chen et al., 2006)

!tag name=Ventilated Chamber Aerosol Deposition (Chen 2006)
     image=../../media/validation/free_flow/isothermal/chen_2006_drift_flux/1_icon.png
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

- A steady, turbulent ventilation flow field in an enclosed chamber using a RANS *k-ε* closure with wall functions.
- Eulerian aerosol transport using a drift–flux formulation that includes:
  - advection by the carrier flow,
  - gravitational settling (drift velocity),
  - Brownian + turbulent diffusion, and
  - near-wall deposition represented as a boundary flux condition.

The test is intended to validate the coupled behavior that matters for indoor/aerosol dispersion problems:
(1) the resolved mean flow field that transports particles through the enclosure and (2) the near-wall mass transfer / deposition behavior that controls concentration gradients and wall losses.

**Out of scope:** thermal buoyancy, humidity effects, coagulation, resuspension, and non-spherical particle effects.

## Problem Description

Chen et al. studied particle distribution and deposition in a mechanically ventilated “scale room” and proposed a drift–flux aerosol model for indoor environments. [!cite](chen2006driftflux)

### Geometry and Operating Conditions (benchmark)

The chamber geometry (Fig. 2 in Chen et al.) is a rectangular enclosure with a square inlet and outlet:

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
$U_{in}=0.225~\mathrm{m/s}$ (reported by Chen et al. as corresponding to 10 air changes per hour). [!cite](chen2006driftflux)

### Aerosol model (benchmark)

Chen et al. model the aerosol phase with a drift–flux transport equation. In the notation used here, a representative form is:

\begin{equation}
\frac{\partial C}{\partial t} + \nabla \cdot \left[(\mathbf{u} + \mathbf{v}_s)\, C \right]
= \nabla \cdot \left(D_\mathrm{eff}\, \nabla C \right) + S_C ,
\end{equation}

where $C$ is particle concentration, $\mathbf{u}$ is the mean air velocity, $\mathbf{v}_s$ is the gravitational settling (drift) velocity,
and $D_\mathrm{eff}$ is the effective diffusivity (Brownian + turbulent/eddy contribution). [!cite](chen2006driftflux)

Wall deposition is represented as a boundary flux using a deposition velocity concept:

\begin{equation}
J_d = v_d\, C_b ,
\end{equation}

where $J_d$ is the deposition flux to the wall, $v_d$ is the deposition velocity, and $C_b$ is a near-wall (boundary layer) concentration. [!cite](chen2006driftflux)

## `OpenPronghorn` Model

This test uses a steady, incompressible finite-volume formulation with a SIMPLE-based segregated pressure–velocity solve and a RANS *k-ε* turbulence model.
Aerosols are modeled as an Eulerian scalar field with drift–flux transport and a wall deposition boundary condition.

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

### Boundary and initial conditions

- **Inlet:** fixed velocity ($U_{in}$) at the inlet patch; turbulence quantities are prescribed from an inlet turbulence intensity and length scale; aerosol concentration is set to $C/C_{in}=1$.
- **Outlet:** outflow condition for velocity and aerosol.
- **Walls:** no-slip with wall functions for turbulence; aerosol wall deposition is imposed via `LinearFVAerosolDepositionBC`.

### Modeled particle properties (this test)

The validation targets the same particle size used for the profile comparisons in Chen et al. (10 µm class). [!cite](chen2006driftflux)

!table id=tab:chen-particles caption=Particle properties used in the `OpenPronghorn` input.
| Property                 |     Symbol |                        Value |
| ------------------------ | ---------: | ---------------------------: |
| Particle diameter        |      $d_p$ | $1.0\times10^{-5}$ m (10 µm) |
| Particle density         |   $\rho_p$ |                1400 kg/m$^3$ |
| Mean free path (air)     |  $\lambda$ |         $6.6\times10^{-8}$ m |
| Turbulent Schmidt number | $Sc_{t,p}$ |                          1.0 |

## Results

The quantities of interest are:

1. **Axial velocity profiles** $u/U_{in}$ vs. height on the chamber center plane at $x=0.2, 0.4, 0.6$ m.
2. **Aerosol concentration profiles** $C/C_{in}$ vs. height at the same locations.

The plots below compare:
- digitized experimental profiles (CSV files in the test `gold/` directory), and
- the corresponding `OpenPronghorn` sampler profiles from the current run.

!media media/validation/aerosols/chen_plot.py
       style=width:95%;margin-left:auto;margin-right:auto;text-align:center
       caption=Axial velocity (top row) and aerosol concentration (bottom row) profiles on the center plane ($y=0.2$ m) at $x=0.2,0.4,0.6$ m. Velocities are normalized by $U_{in}=0.225$ m/s; concentrations are normalized by the inlet concentration $C_{in}$.

## Validation

This test follows the same “error-to-experiment envelope” pattern used by the ERCOFTAC channel-flow validations:

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

## How to run

From the `open_pronghorn` repository root, run the test harness on the aerosol validation directory, for example:

```
./run_tests -j 4 -i validation/aerosols
```

Typical test artifacts include:

- `chen_steady_csv_sampler_line_x_0d*_0002.csv` (current run sampler outputs)
- `gold/` directory containing digitized experiment CSVs and a reference (gold) set of sampler outputs
- `chen_profiles_vnv.py` (validation logic)
- Plots generated by `doc/content/media/validation/aerosols/chen_plot.py` when building documentation.

## Files and directory layout

This test lives under the repository’s aerosol validation area and follows the same structure as other validations (input file, `gold/` reference data, plotting script, and a `validation` configuration file).

## References

!bibtex bibliography=open_pronghorn.bib
