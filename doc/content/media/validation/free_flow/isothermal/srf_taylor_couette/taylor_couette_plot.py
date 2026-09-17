#!/usr/bin/env python3
"""Plot raw Open-Pronghorn Taylor-Couette centroid-profile results.

The input CSV must contain the raw columns written by the centroid-profile
output: ``x``, ``y``, ``vel_abs_x``, ``vel_abs_y``, and ``pressure``.
By default, the repository-relative validation CSV is resolved from the
location of this script. An optional CSV path may still be supplied on the
command line for local testing.
"""

from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Taylor-Couette benchmark parameters
RI = 0.35
RO = 1.0
OMEGA_I = -0.01
OMEGA_O = 0.0
RHO = 1.0
PRESSURE_REFERENCE_RADIUS = 0.75

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CSV_FILE = (
    SCRIPT_DIR / "../../../../../../../validation/free_flow/isothermal/"
    "srf_taylor_couette/"
    "taylor_couette_2d_rel_out_centroid_profile_0001.csv"
).resolve()


def parse_arguments():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv_file",
        nargs="?",
        type=Path,
        default=DEFAULT_CSV_FILE,
        help=(
            "Raw centroid-profile CSV; defaults to the validation-case "
            "output resolved relative to this script"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory for the PNG files (default: script directory)",
    )
    return parser.parse_args()


def resolve_csv(csv_file):
    csv_file = csv_file.expanduser()
    if not csv_file.is_absolute():
        csv_file = (Path.cwd() / csv_file).resolve()
    else:
        csv_file = csv_file.resolve()

    if not csv_file.is_file():
        raise FileNotFoundError(f"CSV file not found: {csv_file}")
    return csv_file


def velocity_coefficients():
    """Return A and B for u_theta(r) = A r + B/r."""
    denominator = RO**2 - RI**2
    coefficient_a = (OMEGA_O * RO**2 - OMEGA_I * RI**2) / denominator
    coefficient_b = (RI**2 * RO**2 * (OMEGA_I - OMEGA_O)) / denominator
    return coefficient_a, coefficient_b


def analytical_velocity(radius):
    """Steady laminar Taylor-Couette tangential velocity."""
    coefficient_a, coefficient_b = velocity_coefficients()
    return coefficient_a * radius + coefficient_b / radius


def analytical_pressure(radius):
    """Taylor-Couette pressure apart from an additive constant."""
    coefficient_a, coefficient_b = velocity_coefficients()
    return RHO * (
        0.5 * coefficient_a**2 * radius**2
        + 2.0 * coefficient_a * coefficient_b * np.log(radius)
        - 0.5 * coefficient_b**2 / radius**2
    )


def load_radial_profile(csv_file):
    """Extract the centroid row closest to the positive x axis."""
    data = pd.read_csv(csv_file)
    required_columns = {
        "x",
        "y",
        "vel_abs_x",
        "vel_abs_y",
        "pressure",
    }
    missing = required_columns.difference(data.columns)
    if missing:
        raise ValueError(
            f"{csv_file.name} is missing required columns: {sorted(missing)}"
        )

    numeric_columns = sorted(required_columns)
    data[numeric_columns] = data[numeric_columns].apply(pd.to_numeric, errors="coerce")
    data = data.dropna(subset=numeric_columns).copy()
    if data.empty:
        raise ValueError(f"{csv_file.name} contains no valid numeric rows")

    data["r"] = np.hypot(data["x"], data["y"])
    data["theta"] = np.arctan2(data["y"], data["x"])

    positive_x = data["x"].to_numpy() > 0.0
    if not np.any(positive_x):
        raise ValueError("No centroid with x > 0 was found in the CSV")

    theta_all = data["theta"].to_numpy()
    candidate_indices = np.flatnonzero(positive_x)
    target_index = candidate_indices[np.argmin(np.abs(theta_all[candidate_indices]))]
    theta_target = theta_all[target_index]

    # Use a wrapped angle difference so this also works near +/-pi.
    angle_difference = np.arctan2(
        np.sin(theta_all - theta_target),
        np.cos(theta_all - theta_target),
    )
    row_mask = np.abs(angle_difference) <= 1.0e-10

    profile = data.loc[row_mask].copy()
    if len(profile) < 2:
        raise ValueError("Could not identify a complete radial centroid row in the CSV")

    # Project Cartesian absolute velocity onto e_theta = (-sin(theta), cos(theta)).
    profile["u_theta"] = (
        -np.sin(profile["theta"]) * profile["vel_abs_x"]
        + np.cos(profile["theta"]) * profile["vel_abs_y"]
    )
    profile = profile.sort_values("r").reset_index(drop=True)
    return profile, theta_target


def build_comparison(profile):
    radius = profile["r"].to_numpy(dtype=float)
    u_theta_numerical = profile["u_theta"].to_numpy(dtype=float)
    pressure_numerical = profile["pressure"].to_numpy(dtype=float)

    if not radius.min() <= PRESSURE_REFERENCE_RADIUS <= radius.max():
        raise ValueError(
            "Pressure reference radius is outside the sampled radial range: "
            f"{radius.min():.6g} <= r <= {radius.max():.6g} m"
        )

    u_theta_analytical = analytical_velocity(radius)
    pressure_analytical = analytical_pressure(radius)

    numerical_reference = np.interp(
        PRESSURE_REFERENCE_RADIUS, radius, pressure_numerical
    )
    analytical_reference = analytical_pressure(PRESSURE_REFERENCE_RADIUS)

    return {
        "r": radius,
        "u_theta_numerical": u_theta_numerical,
        "u_theta_analytical": u_theta_analytical,
        "pressure_numerical_shifted": (pressure_numerical - numerical_reference),
        "pressure_analytical_shifted": (pressure_analytical - analytical_reference),
    }


def plot_velocity(comparison, output_directory):
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(
        comparison["r"],
        comparison["u_theta_numerical"],
        linewidth=2.0,
        label="Open-Pronghorn",
    )
    ax.plot(
        comparison["r"],
        comparison["u_theta_analytical"],
        "--",
        linewidth=2.0,
        label="Analytical",
    )
    ax.set_xlabel("Radius [m]")
    ax.set_ylabel(r"Tangential velocity, $u_\theta$ [m/s]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    output = output_directory / "taylor_couette_velocity_comparison.png"
    fig.savefig(output, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_pressure(comparison, output_directory):
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(
        comparison["r"],
        comparison["pressure_numerical_shifted"],
        linewidth=2.0,
        label="Open-Pronghorn",
    )
    ax.plot(
        comparison["r"],
        comparison["pressure_analytical_shifted"],
        "--",
        linewidth=2.0,
        label="Analytical",
    )
    ax.set_xlabel("Radius [m]")
    ax.set_ylabel("Shifted pressure [Pa]")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()

    output = output_directory / "taylor_couette_pressure_comparison.png"
    fig.savefig(output, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return output


def main():
    arguments = parse_arguments()
    csv_file = resolve_csv(arguments.csv_file)
    output_directory = arguments.output_dir.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    profile, theta_target = load_radial_profile(csv_file)
    comparison = build_comparison(profile)
    velocity_plot = plot_velocity(comparison, output_directory)
    pressure_plot = plot_pressure(comparison, output_directory)

    velocity_error = np.max(
        np.abs(comparison["u_theta_numerical"] - comparison["u_theta_analytical"])
    )
    print(f"Input CSV: {csv_file}")
    print(f"Selected radial row: {np.degrees(theta_target):.8f} deg")
    print(f"Profile points: {len(profile)}")
    print(f"Maximum absolute velocity error: {velocity_error:.6e} m/s")
    print(f"Wrote: {velocity_plot}")
    print(f"Wrote: {pressure_plot}")


if __name__ == "__main__":
    main()
