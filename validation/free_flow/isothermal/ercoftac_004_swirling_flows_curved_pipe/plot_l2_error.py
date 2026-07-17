"""
Compute the combined L2 (and RMS) error in W/W0 (axial velocity) between
OpenPronghorn/MOOSE outputs (v3 and v4) and the ERCOFTAC Anwer & So (case004)
experimental data, at all 9 measurement stations, combining both the
horizontal (beh00, plane AA) and vertical (bev00, plane BB) measurement
planes into one number per station per model version.

Run this AFTER you've executed curvedpipe_noswirl_v3.i and/or
curvedpipe_noswirl_v4.i and have all 18 VectorPostprocessor CSVs (9 stations
x 2 planes) for whichever version(s) you want to check. Missing files are
skipped with a console note rather than causing a crash -- exactly like the
individual per-station plotting scripts.

WHY THE VERTICAL PLANE IS WEIGHTED DOUBLE
------------------------------------------
beh00 (horizontal plane) always measures a FULL diameter, both signs of
r/a, i.e. both the inner and outer bend walls in one profile.
bev00 (vertical plane) only measures a RADIUS (r/a >= 0, the +y half),
relying on the non-swirling case's up/down symmetry to stand in for the
-y half that was never measured.
So a bev00 profile represents twice as much of the pipe cross-section per
data point as a beh00 profile does. To combine the two planes into one
"how wrong is this station overall" number without silently underweighting
the vertical plane, each bev00 SQUARED residual is counted twice before
summing. The one exception is r/a = 0 (the pipe centerline): it is the same
physical point in both planes and has no mirrored partner, so it is
intentionally left un-doubled.

WHY INTERPOLATION, NOT INDEX-MATCHING
--------------------------------------
The simulation's sampler points are placed to match the ERCOFTAC r/a grid
closely, but not necessarily to floating-point-exact agreement, and a
missing/extra near-wall point would silently misalign a simple
index-by-index comparison. Instead, the simulation's W(r/a) curve is
linearly interpolated onto the exact experimental r/a values before
differencing -- this is the standard, robust way to compare a simulation
sampled on one grid against measurements on another.

WHAT "L2 ERROR" MEANS HERE
----------------------------
For a station and a model version:
    L2   = sqrt( sum_i w_i * (W_sim(r_a_i) - W_exp(r_a_i))^2 )
    RMS  = sqrt( sum_i w_i * (...)^2 / sum_i w_i )
where w_i = 1 for every beh00 point and every bev00 point at r/a=0, and
w_i = 2 for every other bev00 point. L2 scales with the (weighted) number
of points combined; RMS does not, so RMS is what's plotted and is directly
comparable in size to ERCOFTAC's own quoted measurement uncertainty on W,
delta(W)/W0 = 0.3%.
"""

import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---- LaTeX-style fonts, matching every other plotting script in this set --
mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
        "axes.grid": False,
        "grid.alpha": 0.3,
        "axes.spines.top": True,
        "axes.spines.right": True,
    }
)

# ---- geometry / normalization constants (must match the .i files) ---------
a = 0.0381          # pipe radius [m]  (D = 0.0762 m)
bulk_u = 10.4       # normalizing bulk velocity used in the .i files [m/s]
Cx, Cz = -0.4953, 0.0   # bend center (x,z)
R = 0.4953              # bend mean radius

# ---- the 9 stations, in physical flow order --------------------------------
# kind='leg'  : centerline_x is the leg's centerline x-coordinate, flow_sign
#               is +1 if flow travels in +z here (inlet leg), -1 if in -z
#               (outlet leg) -- matches every leg script in this set exactly.
# kind='bend' : theta_deg is the bend station angle.
STATIONS = [
    dict(code="sm01",   label=r"$s/D=-1$",          kind="leg",  centerline_x=0.0,      flow_sign=+1),
    dict(code="th0225", label=r"$\theta=22.5^\circ$", kind="bend", theta_deg=22.5),
    dict(code="th0675", label=r"$\theta=67.5^\circ$", kind="bend", theta_deg=67.5),
    dict(code="th1125", label=r"$\theta=112.5^\circ$", kind="bend", theta_deg=112.5),
    dict(code="th1575", label=r"$\theta=157.5^\circ$", kind="bend", theta_deg=157.5),
    dict(code="sp01",   label=r"$s/D=1$",           kind="leg",  centerline_x=-0.9906,  flow_sign=-1),
    dict(code="sp06",   label=r"$s/D=6$",           kind="leg",  centerline_x=-0.9906,  flow_sign=-1),
    dict(code="sp10",   label=r"$s/D=10$",          kind="leg",  centerline_x=-0.9906,  flow_sign=-1),
    dict(code="sp18",   label=r"$s/D=18$",          kind="leg",  centerline_x=-0.9906,  flow_sign=-1),
]

VERSIONS = ["v3", "v4"]
EXP_COLS = ["r_a", "W", "V100", "U100", "wp", "vp", "up", "wu", "uv", "vw"]

# ERCOFTAC's own quoted measurement uncertainty on W (case004 page 4 table),
# used only as a reference line on the plot -- not part of the L2/RMS math.
W_UNCERTAINTY_PCT = 0.3


def load_exp_W(path):
    """Load one ERCOFTAC .dat file's r/a and W columns. Returns None (with a
    console note) if the file doesn't exist, rather than raising."""
    if not path.is_file():
        print(f"    NOTE: experimental file '{path}' not found - skipping.")
        return None
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=EXP_COLS)
    return df[["r_a", "W"]].sort_values("r_a").reset_index(drop=True)


def load_sim_W_beh(path, station):
    """Load one simulation's beh00 (horizontal plane) VPP CSV and compute
    r_a and W/W0 exactly as plot_beh00-*.py does for this station. Returns
    None (with a console note) if the file doesn't exist."""
    if not path.is_file():
        print(f"    NOTE: simulation file '{path}' not found - skipping.")
        return None
    df = pd.read_csv(path)

    if station["kind"] == "leg":
        cx = station["centerline_x"]
        df = df[(df["x"] - cx).abs() < a].copy()
        if station["flow_sign"] > 0:
            df["r_a"] = -df["x"] / a
            df["W"] = df["vel_z"] / bulk_u
        else:
            df["r_a"] = (df["x"] + 0.9906) / a
            df["W"] = -df["vel_z"] / bulk_u
    else:
        theta = np.radians(station["theta_deg"])
        sin_th, cos_th = np.sin(theta), np.cos(theta)
        dist_from_center = ((df["x"] - Cx) ** 2 + (df["z"] - Cz) ** 2) ** 0.5
        df["r_a"] = (R - dist_from_center) / a
        df = df[df["r_a"].abs() < 1.0].copy()
        df["W"] = (-sin_th * df["vel_x"] + cos_th * df["vel_z"]) / bulk_u

    return df[["r_a", "W"]].sort_values("r_a").reset_index(drop=True)


def load_sim_W_bev(path, station):
    """Load one simulation's bev00 (vertical plane) VPP CSV and compute r_a
    and W/W0 exactly as plot_bev00-*.py does for this station. Returns None
    (with a console note) if the file doesn't exist."""
    if not path.is_file():
        print(f"    NOTE: simulation file '{path}' not found - skipping.")
        return None
    df = pd.read_csv(path)
    df = df[df["y"].abs() < a].copy()
    df["r_a"] = df["y"] / a
    df = df[df["r_a"] >= 0].copy()

    if station["kind"] == "leg":
        if station["flow_sign"] > 0:
            df["W"] = df["vel_z"] / bulk_u
        else:
            df["W"] = -df["vel_z"] / bulk_u
    else:
        theta = np.radians(station["theta_deg"])
        sin_th, cos_th = np.sin(theta), np.cos(theta)
        df["W"] = (-sin_th * df["vel_x"] + cos_th * df["vel_z"]) / bulk_u

    return df[["r_a", "W"]].sort_values("r_a").reset_index(drop=True)


def residuals(exp_df, sim_df):
    """Interpolate sim W(r_a) onto the experimental r_a grid and return the
    per-point residual (sim - exp). None in, None out."""
    if exp_df is None or sim_df is None:
        return None
    if len(sim_df) < 2 or len(exp_df) < 1:
        return None
    sim_interp = np.interp(exp_df["r_a"].values, sim_df["r_a"].values, sim_df["W"].values)
    return sim_interp - exp_df["W"].values


def combined_l2(station, version):
    """Combine beh00 + bev00 residuals for one station/version into a single
    weighted L2 norm and RMS, per the weighting rule in the module docstring."""
    beh_exp = load_exp_W(Path(f"swb-allfiles/beh00-{station['code']}.dat"))
    bev_exp = load_exp_W(Path(f"swb-allfiles/bev00-{station['code']}.dat"))
    beh_sim = load_sim_W_beh(Path(f"curvedpipe_noswirl_{version}_csv_beh00_{station['code']}_0002.csv"), station)
    bev_sim = load_sim_W_bev(Path(f"curvedpipe_noswirl_{version}_csv_bev00_{station['code']}_0002.csv"), station)

    beh_res = residuals(beh_exp, beh_sim)
    bev_res = residuals(bev_exp, bev_sim)

    if beh_res is None and bev_res is None:
        return None

    weighted_sq_sum = 0.0
    n_eff = 0.0

    if beh_res is not None:
        weighted_sq_sum += np.sum(beh_res ** 2)
        n_eff += len(beh_res)

    if bev_res is not None:
        # Double every bev00 point EXCEPT r/a=0, which is the same physical
        # point as beh00's own r/a=0 and has no mirrored -y partner.
        r_a_bev = bev_exp["r_a"].values
        weights = np.where(r_a_bev == 0.0, 1.0, 2.0)
        weighted_sq_sum += np.sum(weights * bev_res ** 2)
        n_eff += np.sum(weights)

    l2 = np.sqrt(weighted_sq_sum)
    rms = np.sqrt(weighted_sq_sum / n_eff) if n_eff > 0 else np.nan
    return dict(l2=l2, rms=rms, n_eff=n_eff,
                n_beh=(len(beh_res) if beh_res is not None else 0),
                n_bev=(len(bev_res) if bev_res is not None else 0))


# ---- 1. compute everything ------------------------------------------------
print("Computing W/W0 L2 error by station (beh00 + weighted bev00)...\n")
results = {}
for station in STATIONS:
    # print(f"  {station['code']}")
    results[station["code"]] = {}
    for version in VERSIONS:
        results[station["code"]][version] = combined_l2(station, version)

# ---- 2. print console table -------------------------------------------------
rows = []
for station in STATIONS:
    row = {"Station": station["label"].replace(r"\theta", "theta").replace("$", "").replace(r"\circ", "deg")}
    for version in VERSIONS:
        r = results[station["code"]][version]
        row[f"{version} L2"] = r["l2"] if r is not None else np.nan
        row[f"{version} RMS %"] = 100.0 * r["rms"] if r is not None else np.nan
    r_any = results[station["code"]]["v3"] or results[station["code"]]["v4"]
    row["N"] = r_any["n_eff"] if r_any is not None else np.nan
    rows.append(row)

table = pd.DataFrame(rows).set_index("Station")
pd.set_option("display.float_format", lambda x: f"{x:,.4f}")
print("\n" + "=" * 78)
print("W/W0 combined L2 error (beh00 + 2x-weighted bev00) vs. ERCOFTAC")
print("=" * 78)
print(table.to_string())
print("=" * 78)

# ---- 3. bar chart -----------------------------------------------------------
codes = [s["code"] for s in STATIONS]
labels = [s["label"] for s in STATIONS]
rms_v3 = [100.0 * results[c]["v3"]["rms"] if results[c]["v3"] is not None else np.nan for c in codes]
rms_v4 = [100.0 * results[c]["v4"]["rms"] if results[c]["v4"] is not None else np.nan for c in codes]

x = np.arange(len(codes))
width = 0.35

fig, ax = plt.subplots(figsize=(12, 5.5))

bars_v3 = ax.bar(x - width / 2, rms_v3, width, color="#1f77b4", label="MOOSE v3")
bars_v4 = ax.bar(x + width / 2, rms_v4, width, color="#d62728", label="MOOSE v4")

ax.axhline(W_UNCERTAINTY_PCT, color="black", linestyle="--", linewidth=0.9,
           label=f"ERCOFTAC measurement uncertainty ({W_UNCERTAINTY_PCT}%)")

for bars in (bars_v3, bars_v4):
    for b in bars:
        h = b.get_height()
        if np.isfinite(h):
            ax.annotate(f"{h:.2f}", xy=(b.get_x() + b.get_width() / 2, h),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8)

ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel(r"RMS error in $W/W_0$  (%)")
ax.set_xlabel("Station")
ax.set_title(r"Combined $W/W_0$ error vs. ERCOFTAC, by station")
ax.grid(alpha=0.3, axis="y")
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.0)
# headroom above the tallest bar so the value labels never get clipped, and
# the legend sits outside the axes so it never overlaps a bar or its label.
ax.set_ylim(0, max(np.nanmax(rms_v3 + rms_v4), 1.0) * 1.18)
# data series first, reference line last, regardless of draw order above.
handles, leg_labels = ax.get_legend_handles_labels()
order = [1, 2, 0]
ax.legend([handles[i] for i in order], [leg_labels[i] for i in order],
          frameon=True, fontsize=9, handlelength=2.0,
          loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)

fig.tight_layout()
fig.savefig("l2_error_W.png", dpi=300, bbox_inches="tight")
print("\nSaved l2_error_W.png")