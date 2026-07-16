"""
Compare OpenPronghorn/MOOSE outputs (v3 and v4 input files) against the
ERCOFTAC Anwer & So (case004) experimental data, at s/D = -1 upstream of
the bend, vertical plane (plane BB). Only the +y half of each simulation is
plotted (see the geometry note in the .i files: y is the radial coordinate).

Run this AFTER you've executed curvedpipe_noswirl_v3.i and/or
curvedpipe_noswirl_v4.i and have their VectorPostprocessor CSV output in
hand. Either simulation CSV may be absent - this script will flag that in
the console and just skip it, rather than fail. The experimental .dat file
is required, though: if it's missing, this script stops with an error,
since there's nothing meaningful to plot without it.
"""

import sys
from pathlib import Path

import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt

# ---- LaTeX-style fonts, matching 2d_channel_forch_plot_line_samples.py ----
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

# ---- 1. load experimental data (required - fail loudly if missing) --------
exp_path = Path('swb-allfiles/bev00-sm01.dat')
if not exp_path.is_file():
    sys.exit(f"ERROR: experimental data file not found: {exp_path}\n"
             f"Nothing to plot without it - stopping.")

exp_cols = ['r_a', 'W', 'V100', 'U100', 'wp', 'vp', 'up', 'wu', 'uv', 'vw']
exp = pd.read_csv(exp_path, sep=r'\s+', comment='#', header=None, names=exp_cols)
exp['V'] = exp['V100'] / 100.0   # file stores V, U scaled by 100
exp['U'] = exp['U100'] / 100.0

# ---- 2. load simulation data (optional - flag and skip if missing) --------
v4_csv = 'curvedpipe_noswirl_v4_csv_bev00_sm01_0002.csv'
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

    # Exclude the two endpoint samples (r/a = +-1, i.e. y = +-a exactly).
    # These aren't real near-wall data: on a wall-function mesh the first
    # cell center sits away from the wall, so the sampler returns that
    # cell's value at the wall point instead of the enforced (zero) no-slip
    # value. The experiment sidesteps this by never measuring past
    # r/a = 0.938, so we drop the endpoints here too.
    sim = sim[sim['y'].abs() < a].copy()

    # The sampled line runs along y through the pipe center, so y IS the
    # radial coordinate. Keep only the +y half, per this round's request.
    sim['r_a'] = sim['y'] / a
    sim = sim[sim['r_a'] >= 0].sort_values('r_a')

    # Flow in this leg travels in +z here (this is the INLET leg: flow enters
    # at the inlet boundary z=-0.1524 and moves toward the bend entrance at
    # z=0 -- confirmed directly from the mesh file's inlet/outlet side-set
    # node coordinates, not just the .i file's comment). An earlier version
    # of this comment said "-z", which was wrong; the code below was already
    # correct (unnegated vel_z), so no numerical change here, comment only.
    sim['W'] = sim['vel_z'] / bulk_u
    sim['U'] = sim['vel_y'] / bulk_u   # radial component, along the sampled line
    # CORRECTED SIGN: proper right-handed (radial, circumferential, axial)
    # cylindrical convention requires V_hat = W_hat x U_hat, not U_hat x W_hat.
    # At this leg (W=+z_hat, U=+y_hat): W x U = z_hat x y_hat = -x_hat.
    # So positive V = toward the INNER bend wall here (an earlier version of
    # this line used +vel_x, which was backwards).
    sim['V'] = -sim['vel_x'] / bulk_u   # circumferential component, in-plane; positive = toward inner wall
    return sim


v3 = load_sim(v3_csv)
v4 = load_sim(v4_csv)

# ---- published measurement uncertainty (ERCOFTAC case004, page 4 table) ---
# Absolute error bands, in the same W0-normalized units as the data itself -
# not a percentage of the local value (see delta(W)/Wo etc. in the source doc).
uncertainty = {'W': 0.003, 'U': 0.03, 'V': 0.07}

# ---- 3. plot ----------------------------------------------------------------
# Points only, no connecting lines between samples - the format strings below
# deliberately carry no linestyle character (matplotlib defaults to a marker
# with no connecting line when none is given).
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
fig.suptitle(r's/D = -1, vertical plane (plane BB)')
fig.tight_layout()
fig.savefig('comparison_bev00-sm01.png', dpi=300, bbox_inches='tight')
print('Saved comparison_bev00-sm01.png')