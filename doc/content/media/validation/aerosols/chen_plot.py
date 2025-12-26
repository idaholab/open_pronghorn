import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

"""
Publication-ready Chen validation plot generator (MOOSE-friendly).

Key requirement from your "old script":
- Velocity: plot experimental as-is (already nondimensional), and plot simulation as vel_x / U_REF
- y-axis: z in meters, y-limits [0, H]
- Concentration: plot exp as-is, sim as aerosol / C_REF, with x-limits [0, 1]

This script reads everything from a validation "gold" directory.
It saves a single PNG next to this script (default: chen_profiles.png).
"""

# Geometry
H = 0.4  # m

# Nondimensional references (MATCH OLD SCRIPT BEHAVIOR)
U_REF = 1.0
C_REF = 1.0


def _set_publication_style():
    plt.rcParams.update(
        {
            "savefig.dpi": 600,
            "font.size": 10.5,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def _apply_axes_style(ax):
    ax.grid(True, which="major", linewidth=0.6, alpha=0.25)
    ax.grid(True, which="minor", linewidth=0.4, alpha=0.15)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))


def load_exp_two_column_csv(path: Path):
    """
    Old-script behavior:
      - First column = non-dimensional dependent variable (u/Uref or C/Cref)
      - Second column = z in meters
    """
    df = pd.read_csv(path)

    # Allow placeholder/empty files during early setup
    if df.shape[0] == 0:
        return np.array([]), np.array([])

    dep = df[df.columns[0]].to_numpy(dtype=float)
    z = df[df.columns[1]].to_numpy(dtype=float)
    return dep, z


def load_sim_profile_from_sampler(path: Path, value_name: str, ref_val: float):
    """
    Old-script behavior:
      - z is taken as-is (meters) from sampler
      - value is divided by ref_val
    """
    df = pd.read_csv(path, comment="#")
    if "z" not in df.columns:
        raise RuntimeError(f"Could not find 'z' column in {path}")
    if value_name not in df.columns:
        raise RuntimeError(f"Could not find '{value_name}' column in {path}")

    z = df["z"].to_numpy(dtype=float)
    val = df[value_name].to_numpy(dtype=float)

    idx = np.argsort(z)
    z = z[idx]
    val = val[idx]

    return val / ref_val, z


def add_exp_with_errorbars(ax, dep, z, label, color, marker, fixed=False):
    dep = np.asarray(dep, dtype=float)
    z = np.asarray(z, dtype=float)
    if dep.size == 0:
        return

    xerr = 0.1 * np.max(np.abs(dep))
    if fixed:
        xerr = 0.08

    ax.errorbar(
        dep,
        z,
        xerr=xerr,
        yerr=None,
        label=label,
        marker=marker,
        linestyle="none",
        markersize=5.0,
        markeredgewidth=0.9,
        color=color,
        ecolor=color,
        elinewidth=0.9,
        capsize=2.0,
        capthick=0.9,
        zorder=3,
    )


def add_sim_line(ax, dep_nd, z_m, label, color, linestyle="-"):
    ax.plot(dep_nd, z_m, linestyle=linestyle, marker='o', color=color, linewidth=0.5, label=label, zorder=2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--valdir",
        default=None,
        help="Directory containing the gold/ CSV files. If omitted, uses ../../../../../validation/aerosols/gold relative to this script.",
    )
    parser.add_argument(
        "--outfile",
        default="chen_profiles.png",
        help="Output PNG filename (saved next to this script unless absolute).",
    )
    parser.add_argument("--dpi", type=int, default=600, help="PNG DPI (default: 600).")
    parser.add_argument("--no-show", action="store_true", help="Do not call plt.show() (useful for MOOSE batch).")
    args = parser.parse_args()

    _set_publication_style()

    script_dir = Path(__file__).resolve().parent

    # Default to the same type of relative path you used in your working script
    if args.valdir is None:
        val_dir = (script_dir / "../../../../../validation/aerosols/gold").resolve()
    else:
        val_dir = Path(args.valdir).expanduser().resolve()

    out_path = Path(args.outfile)
    if not out_path.is_absolute():
        out_path = (script_dir / out_path).resolve()

    # ---- Experimental (gold) ----
    u_exp_0d2, z_u_0d2 = load_exp_two_column_csv(val_dir / "data_u_x_0d2.csv")
    u_exp_0d4, z_u_0d4 = load_exp_two_column_csv(val_dir / "data_u_x_0d4.csv")
    u_exp_0d6, z_u_0d6 = load_exp_two_column_csv(val_dir / "data_u_x_0d6.csv")

    c_exp_0d2, z_c_0d2 = load_exp_two_column_csv(val_dir / "data_c_x_0d2.csv")
    c_exp_0d4, z_c_0d4 = load_exp_two_column_csv(val_dir / "data_c_x_0d4.csv")
    c_exp_0d6, z_c_0d6 = load_exp_two_column_csv(val_dir / "data_c_x_0d6.csv")

    # ---- Simulation sampler output (gold or current-run location you choose) ----
    u_sim_0d2, z_sim_0d2 = load_sim_profile_from_sampler(val_dir / "chen_steady_csv_sampler_line_x_0d2_0002.csv", "vel_x", U_REF)
    u_sim_0d4, z_sim_0d4 = load_sim_profile_from_sampler(val_dir / "chen_steady_csv_sampler_line_x_0d4_0002.csv", "vel_x", U_REF)
    u_sim_0d6, z_sim_0d6 = load_sim_profile_from_sampler(val_dir / "chen_steady_csv_sampler_line_x_0d6_0002.csv", "vel_x", U_REF)

    c_sim_0d2, zc_sim_0d2 = load_sim_profile_from_sampler(val_dir / "chen_steady_csv_sampler_line_x_0d2_0002.csv", "aerosol", C_REF)
    c_sim_0d4, zc_sim_0d4 = load_sim_profile_from_sampler(val_dir / "chen_steady_csv_sampler_line_x_0d4_0002.csv", "aerosol", C_REF)
    c_sim_0d6, zc_sim_0d6 = load_sim_profile_from_sampler(val_dir / "chen_steady_csv_sampler_line_x_0d6_0002.csv", "aerosol", C_REF)

    # ---- Figure ----
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8), sharey=True)
    (ax_u_0d2, ax_u_0d4, ax_u_0d6), (ax_c_0d2, ax_c_0d4, ax_c_0d6) = axes

    # Panel labels
    panel_labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    for ax, lab in zip(axes.ravel(), panel_labels):
        ax.text(0.02, 0.98, lab, transform=ax.transAxes, va="top", ha="left",
                fontsize=11, fontweight="bold")

    # ---- Velocity row (LIKE YOUR OLD SCRIPT) ----
    add_exp_with_errorbars(ax_u_0d2, u_exp_0d2, z_u_0d2, "Exp", "black", "o")
    add_sim_line(ax_u_0d2, u_sim_0d2, z_sim_0d2, "Sim", "tab:blue")
    ax_u_0d2.set_xlabel(r"$u / U_{\mathrm{ref}}$")
    ax_u_0d2.set_ylabel(r"$z$ [m]")
    ax_u_0d2.set_title(r"Axial velocity, $x = 0.2$ m", loc="center")

    add_exp_with_errorbars(ax_u_0d4, u_exp_0d4, z_u_0d4, "Exp", "black", "o")
    add_sim_line(ax_u_0d4, u_sim_0d4, z_sim_0d4, "Sim", "tab:blue")
    ax_u_0d4.set_xlabel(r"$u / U_{\mathrm{ref}}$")
    ax_u_0d4.set_title(r"Axial velocity, $x = 0.4$ m", loc="center")

    add_exp_with_errorbars(ax_u_0d6, u_exp_0d6, z_u_0d6, "Exp", "black", "s")
    add_sim_line(ax_u_0d6, u_sim_0d6, z_sim_0d6, "Sim", "tab:blue")
    ax_u_0d6.set_xlabel(r"$u / U_{\mathrm{ref}}$")
    ax_u_0d6.set_title(r"Axial velocity, $x = 0.6$ m", loc="center")

    # ---- Concentration row (LIKE YOUR OLD SCRIPT) ----
    add_exp_with_errorbars(ax_c_0d2, c_exp_0d2, z_c_0d2, "Exp", "black", "^", fixed=True)
    add_sim_line(ax_c_0d2, c_sim_0d2, zc_sim_0d2, "Sim", "tab:red")
    ax_c_0d2.set_xlabel(r"$C / C_{\mathrm{ref}}$")
    ax_c_0d2.set_ylabel(r"$z$ [m]")
    ax_c_0d2.set_title(r"Aerosol conc., $x = 0.2$ m", loc="center")
    ax_c_0d2.set_xlim(0.0, 1.0)

    add_exp_with_errorbars(ax_c_0d4, c_exp_0d4, z_c_0d4, "Exp", "black", "^", fixed=True)
    add_sim_line(ax_c_0d4, c_sim_0d4, zc_sim_0d4, "Sim", "tab:red")
    ax_c_0d4.set_xlabel(r"$C / C_{\mathrm{ref}}$")
    ax_c_0d4.set_title(r"Aerosol conc., $x = 0.4$ m", loc="center")
    ax_c_0d4.set_xlim(0.0, 1.0)

    add_exp_with_errorbars(ax_c_0d6, c_exp_0d6, z_c_0d6, "Exp", "black", "v", fixed=True)
    add_sim_line(ax_c_0d6, c_sim_0d6, zc_sim_0d6, "Sim", "tab:red")
    ax_c_0d6.set_xlabel(r"$C / C_{\mathrm{ref}}$")
    ax_c_0d6.set_title(r"Aerosol conc., $x = 0.6$ m", loc="center")
    ax_c_0d6.set_xlim(0.0, 1.0)

    # Common styling
    for ax in axes.ravel():
        ax.set_ylim(0.0, H)
        _apply_axes_style(ax)

    # Legends: only left-most in each row to reduce clutter
    ax_u_0d2.legend(frameon=False, loc="lower right")
    ax_c_0d2.legend(frameon=False, loc="lower right")

    fig.tight_layout(pad=0.9)

    fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote: {out_path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
