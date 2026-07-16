"""
Compare OpenPronghorn/MOOSE outputs (v3 and v4 input files) against the
ERCOFTAC Anwer & So (case004) experimental data, at theta = 22.5 degrees
within the bend, vertical plane (plane BB). Only the +y half of each
simulation is plotted (matches the .dat file, which is single-sided here too).

Run this AFTER you've executed curvedpipe_noswirl_v3.i and/or
curvedpipe_noswirl_v4.i and have their VectorPostprocessor CSV output in
hand. Either simulation CSV may be absent - this script will flag that in
the console and just skip it, rather than fail. The experimental .dat file
is required, though: if it's missing, this script stops with an error,
since there's nothing meaningful to plot without it.

GEOMETRY NOTE (new relative to the sm01/sp01 scripts): at a bend station the
local flow-aligned frame is ROTATED relative to global (x,y,z). The axial
("W") direction is the local tangent T(theta) = (-sin theta, 0, cos theta),
not simply vel_z. y remains the vertical/out-of-plane direction unchanged
(no rotation needed there), matching the sm01/sp01 convention exactly at the
theta=0/180 limits:
    theta=0   -> W = vel_z         (matches plot_bev00-sm01.py exactly)
    theta=180 -> W = -vel_z        (matches plot_bev00-sp01.py exactly)

CAVEAT on the "V" (secondary, in-plane) component: the sm01/sp01 scripts you
already have both use V = +vel_x with NO rotation and no sign flip between
legs (this was flagged earlier in our conversation as an unresolved
simplification, not a verified convention). The geometrically consistent
extension to a rotated bend station is V = cos(theta)*vel_x + sin(theta)*vel_z
(the raw in-plane-normal projection), which correctly reduces to +vel_x at
theta=0 (matching sm01) but reduces to -vel_x at theta=180 -- the OPPOSITE
sign from what the given sp01 script uses. That's not a bug here; it's the
existing sp01 script's V that isn't rotated. I've implemented the physically
consistent rotation below and flagged it -- let me know if you'd rather I
revisit the leg scripts to match instead.
"""

import sys
from pathlib import Path

import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt

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
a = 0.0381        # pipe radius [m]  (D = 0.0762 m)
bulk_u = 10.4     # normalizing bulk velocity used in the .i files [m/s]

# ---- bend station angle (verified against wall mesh in an earlier turn) ---
theta_deg = 22.5
sin_th = 0.3826834324
cos_th = 0.9238795325

# ---- 1. load experimental data (required - fail loudly if missing) --------
exp_path = Path('swb-allfiles/bev00-th0225.dat')
if not exp_path.is_file():
    sys.exit(f"ERROR: experimental data file not found: {exp_path}\n"
             f"Nothing to plot without it - stopping.")

exp_cols = ['r_a', 'W', 'V100', 'U100', 'wp', 'vp', 'up', 'wu', 'uv', 'vw']
exp = pd.read_csv(exp_path, sep=r'\s+', comment='#', header=None, names=exp_cols)
exp['V'] = exp['V100'] / 100.0   # file stores V, U scaled by 100
exp['U'] = exp['U100'] / 100.0

# ---- 2. load simulation data (optional - flag and skip if missing) --------
v4_csv = 'curvedpipe_noswirl_v4_csv_bev00_th0225_0002.csv'
v3_csv = v4_csv.replace('_v4_', '_v3_')   # assumed to be named the same way as v4


def load_sim(csv_path):
    """Load one simulation's VectorPostprocessor CSV and convert it into the
    experiment's coordinates/units. Returns None (and prints a note) if the
    file doesn't exist, rather than raising."""
    path = Path(csv_path)
    if not path.is_file():
        print(f"NOTE: '{csv_path}' not found - skipping this series.")
        return None

    sim = pd.read_csv(path)

    # Exclude the two endpoint samples (r/a = +-1, i.e. y = +-a exactly) --
    # same reasoning as the sm01/sp01 scripts.
    sim = sim[sim['y'].abs() < a].copy()

    # y is still the radial coordinate here (unrotated) -- keep +y half only.
    sim['r_a'] = sim['y'] / a
    sim = sim[sim['r_a'] >= 0].sort_values('r_a')

    # Rotated axial direction: W = velocity dot T(theta), T(theta) =
    # (-sin theta, 0, cos theta). Reduces exactly to sm01's vel_z at
    # theta=0 and to sp01's -vel_z at theta=180.
    sim['W'] = (-sin_th * sim['vel_x'] + cos_th * sim['vel_z']) / bulk_u
    sim['U'] = sim['vel_y'] / bulk_u   # vertical direction never rotates
    # Rotated in-plane ("circumferential") component -- see the CAVEAT in
    # the module docstring above regarding the sign mismatch with sp01/sm01
    # at the theta=180 limit.
    sim['V'] = (cos_th * sim['vel_x'] + sin_th * sim['vel_z']) / bulk_u
    return sim


v3 = load_sim(v3_csv)
v4 = load_sim(v4_csv)

# ---- published measurement uncertainty (ERCOFTAC case004, page 4 table) ---
uncertainty = {'W': 0.003, 'U': 0.03, 'V': 0.07}

# ---- 3. plot ----------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, quantity, ylabel in zip(
    axes,
    ['W', 'U', 'V'],
    [r'$W / W_0$  (axial)', r'$U / W_0$  (radial)', r'$V / W_0$  (circumferential)'],
):
    ax.errorbar(exp['r_a'], exp[quantity], yerr=uncertainty[quantity],
                fmt='o', color='black', markersize=4.0, capsize=3,
                elinewidth=0.8, capthick=0.8, label='Experiment')
    if v3 is not None:
        ax.plot(v3['r_a'], v3[quantity], 's', color='#1f77b4', markersize=3.2, label='MOOSE v3')
    if v4 is not None:
        ax.plot(v4['r_a'], v4[quantity], '^', color='#d62728', markersize=3.2, label='MOOSE v4')
    ax.set_xlabel(r'$r / a$')
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)

axes[0].legend(frameon=True, fontsize=8, handlelength=2.0)
fig.suptitle(r'theta = 22.5 deg, vertical plane (plane BB)')
fig.tight_layout()
fig.savefig('comparison_bev00-th0225.png', dpi=300, bbox_inches='tight')
print('Saved comparison_bev00-th0225.png')
