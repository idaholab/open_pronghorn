"""
Compare OpenPronghorn/MOOSE outputs (v3 and v4 input files) against the
ERCOFTAC Anwer & So (case004) experimental data, at theta = 22.5 degrees
within the bend, horizontal plane (plane AA) -- the FULL DIAMETER (both
inner and outer bend sides), same as the other beh00 scripts.

Run this AFTER you've executed curvedpipe_noswirl_v3.i and/or
curvedpipe_noswirl_v4.i and have their VectorPostprocessor CSV output in
hand. Either simulation CSV may be absent - this script will flag that in
the console and just skip it, rather than fail. The experimental .dat file
is required, though: if it's missing, this script stops with an error,
since there's nothing meaningful to plot without it.

GEOMETRY: at a bend station, recovering r/a from the sampled (x,y,z) uses a
theta-INDEPENDENT formula -- r_a = (R - dist_from_bend_center) / a, where
dist_from_bend_center = sqrt((x-Cx)**2 + (z-Cz)**2). This only holds on the
curved section itself (not the straight legs, which use their own simpler
formula in the other beh00 scripts), but it doesn't need theta explicitly,
which makes it robust to any small numerical noise in the sampled point
positions. Sign matches the same Dean-effect convention validated for the
leg stations: negative r/a = outer wall (high W), positive r/a = inner wall.

W and U (the along-sweep, in-bend-plane radial component) both need the same
theta-rotation as the vertical bend scripts:
    W = velocity . T(theta),  T(theta) = (-sin theta, 0, cos theta)
    U = -(velocity . N(theta)),  N(theta) = (cos theta, 0, sin theta)
(the minus sign on U matches the same "positive U = toward inner wall,
increasing r_a" convention used in the leg horizontal scripts -- this is a
direct, derived consequence of the r_a sign convention, not a new guess).
V (vertical, out-of-plane) is simply vel_y, unchanged, same as every other
beh00 script -- the vertical direction never rotates.
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

# ---- bend geometry (verified against the wall mesh in an earlier turn) ----
Cx, Cz = -0.4953, 0.0
R = 0.4953
theta_deg = 22.5
sin_th = 0.3826834324
cos_th = 0.9238795325

# ---- 1. load experimental data (required - fail loudly if missing) --------
exp_path = Path('swb-allfiles/beh00-th0225.dat')
if not exp_path.is_file():
    sys.exit(f"ERROR: experimental data file not found: {exp_path}\n"
             f"Nothing to plot without it - stopping.")

exp_cols = ['r_a', 'W', 'V100', 'U100', 'wp', 'vp', 'up', 'wu', 'uv', 'vw']
exp = pd.read_csv(exp_path, sep=r'\s+', comment='#', header=None, names=exp_cols)
exp['V'] = exp['V100'] / 100.0   # file stores V, U scaled by 100
exp['U'] = exp['U100'] / 100.0

# ---- 2. load simulation data (optional - flag and skip if missing) --------
v4_csv = 'curvedpipe_noswirl_v4_csv_beh00_th0225_0002.csv'
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

    # Distance from the bend center, and the theta-independent r_a formula.
    dist_from_center = ((sim['x'] - Cx)**2 + (sim['z'] - Cz)**2)**0.5
    sim['r_a'] = (R - dist_from_center) / a

    # Exclude the two endpoint samples (r/a = +-1 exactly) -- same reasoning
    # as every other script here.
    sim = sim[sim['r_a'].abs() < 1.0].copy()
    sim = sim.sort_values('r_a')

    # Rotated axial (W) and in-plane-radial (U) components -- see the module
    # docstring above for the derivation.
    sim['W'] = (-sin_th * sim['vel_x'] + cos_th * sim['vel_z']) / bulk_u
    sim['U'] = (cos_th * sim['vel_x'] + sin_th * sim['vel_z']) / bulk_u   # positive = toward outer bend wall
    sim['V'] = sim['vel_y'] / bulk_u   # vertical direction, never rotates
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
    [r'$W / W_0$  (axial)', r'$U / W_0$  (radial)', r'$V / W_0$  (vertical)'],
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
fig.suptitle(r'theta = 22.5 deg, horizontal plane (plane AA)')
fig.tight_layout()
fig.savefig('comparison_beh00-th0225.png', dpi=300, bbox_inches='tight')
print('Saved comparison_beh00-th0225.png')
