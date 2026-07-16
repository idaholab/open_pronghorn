"""
Compare OpenPronghorn/MOOSE outputs (v3 and v4 input files) against the
ERCOFTAC Anwer & So (case004) experimental data, at s/D = -1 upstream of the bend,
horizontal plane (plane AA) -- the FULL DIAMETER, both inner and outer bend
sides, unlike the vertical (bev00) scripts which only use the +y half.

Run this AFTER you've executed curvedpipe_noswirl_v3.i and/or
curvedpipe_noswirl_v4.i and have their VectorPostprocessor CSV output in
hand. Either simulation CSV may be absent - this script will flag that in
the console and just skip it, rather than fail. The experimental .dat file
is required, though: if it's missing, this script stops with an error,
since there's nothing meaningful to plot without it.

SIGN CONVENTION (verified against the Dean-effect W asymmetry in your own
beh00-sp01.dat / beh00-th1575.dat -- negative r/a is the OUTER bend wall,
where axial velocity is elevated; positive r/a is the INNER wall):
    inlet leg  (x=0 centerline):       r_a = -x / a
    outlet leg (x=-0.9906 centerline): r_a = (x + 0.9906) / a
The corresponding "U" (along-the-sweep, i.e. radial-in-bend-plane) velocity
component is defined with the SAME sign convention as r_a itself (positive U
= flow toward increasing r_a = toward the inner wall), which is why its sign
differs between the inlet and outlet legs below -- this isn't an assumption,
it falls directly out of differentiating the r_a formula above.
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
centerline_x = 0.0   # this leg's centerline x-coordinate

# ---- 1. load experimental data (required - fail loudly if missing) --------
exp_path = Path('swb-allfiles/beh00-sm01.dat')
if not exp_path.is_file():
    sys.exit(f"ERROR: experimental data file not found: {exp_path}\n"
             f"Nothing to plot without it - stopping.")

exp_cols = ['r_a', 'W', 'V100', 'U100', 'wp', 'vp', 'up', 'wu', 'uv', 'vw']
exp = pd.read_csv(exp_path, sep=r'\s+', comment='#', header=None, names=exp_cols)
exp['V'] = exp['V100'] / 100.0   # file stores V, U scaled by 100
exp['U'] = exp['U100'] / 100.0

# ---- 2. load simulation data (optional - flag and skip if missing) --------
v4_csv = 'curvedpipe_noswirl_v4_csv_beh00_sm01_0002.csv'
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

    # Exclude the two endpoint samples (r/a = +-1 exactly) -- same reasoning
    # as the vertical scripts: the wall-adjacent cell center isn't the true
    # no-slip wall value, and the experiment never measures past r/a=0.938.
    sim = sim[(sim['x'] - centerline_x).abs() < a].copy()

    # This is a FULL DIAMETER station (both signs of r_a) -- no +side-only
    # filtering here, unlike the vertical scripts.
    sim['r_a'] = -sim['x'] / a
    sim = sim.sort_values('r_a')

    # Flow in this leg travels in +z (per the .i file's geometry comment).
    sim['W'] = sim['vel_z'] / bulk_u
    sim['U'] = sim['vel_x'] / bulk_u   # radial (in-bend-plane); positive = toward outer bend wall
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
fig.suptitle(r's/D = -1, horizontal plane (plane AA)')
fig.tight_layout()
fig.savefig('comparison_beh00-sm01.png', dpi=300, bbox_inches='tight')
print('Saved comparison_beh00-sm01.png')
