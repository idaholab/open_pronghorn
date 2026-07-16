#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import pandas as pd

os.chdir(os.path.dirname(os.path.realpath(__file__)))

csv_file = "../../../../../../../../validation/free_flow/isothermal/vortex_shedding/cylinder/gold/force_coefficients.csv"
data = pd.read_csv(csv_file).iloc[0]

quantities = (
    {
        "column": "Maximum drag coefficient",
        "label": r"Maximum drag coefficient, $C_{D,\max}$",
        "bounds": (3.22, 3.24),
    },
    {
        "column": "Maximum lift coefficient",
        "label": r"Maximum lift coefficient, $C_{L,\max}$",
        "bounds": (0.99, 1.01),
    },
)

plt.rcParams.update({"font.family": "serif", "font.size": 12})
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

for ax, quantity in zip(axes, quantities):
    value = data[quantity["column"]]
    lower_bound, upper_bound = quantity["bounds"]
    bound_width = upper_bound - lower_bound

    ax.axhspan(
        lower_bound,
        upper_bound,
        color="#56B4E9",
        alpha=0.3,
        label="Acceptable range",
    )
    ax.plot(1, value, "X", color="#D55E00", markersize=10, label="Current value")
    ax.set_xlim(0.5, 1.5)
    ax.set_ylim(
        min(value, lower_bound) - 0.5 * bound_width,
        max(value, upper_bound) + 0.5 * bound_width,
    )
    ax.get_xaxis().set_visible(False)
    ax.set_ylabel(quantity["label"])
    ax.grid(True, which="major", linestyle="-", linewidth=0.5, color="grey")
    ax.legend()

fig.tight_layout()
plt.savefig("force_coefficients.png", dpi=300)
plt.close()
