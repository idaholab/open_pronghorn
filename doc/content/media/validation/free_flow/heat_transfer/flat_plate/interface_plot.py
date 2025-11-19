import matplotlib.pyplot as plt
import os
import pandas as pd
from common import (
    read_moose,
    analytic_approximate_solution_dht_erfcx,
    analytic_approximate_solution_bl,
)

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Formal appearance
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 12,
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
    }
)

# Read the files
print_commercial_cfd = False
moose_fluid_df = read_moose(
    "../../../../../../../validation/free_flow/heat_transfer/flat_plate/gold/cht_rob-rob_out_interface_temp_0001.csv",
    "x",
)
moose_fluid_df = moose_fluid_df.sort_values("x")
if print_commercial_cfd:
    reference_cfd_df = read_moose(
        "../../../../../../../validation/free_flow/heat_transfer/flat_plate/gold/reference_cfd_horizontal.csv"
    ).sort_values("x")

# Compute analytic results
st_bl = [analytic_approximate_solution_bl(x, 1e-9) for x in moose_fluid_df["x"]]
st_dht = [analytic_approximate_solution_dht_erfcx(x, 1e-9) for x in moose_fluid_df["x"]]

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(moose_fluid_df["x"], moose_fluid_df["T"], "b", label="MOOSE/Pronghorn")
ax.plot(moose_fluid_df["x"], st_bl, "r", label="Analytic BL (Luikov)")
ax.plot(moose_fluid_df["x"], st_dht, "k", label="Analytic DHT (Luikov)")

if print_commercial_cfd:
    ax.plot(
        reference_cfd_df["x"], reference_cfd_df["T"], "--", label="Commercial CFD tool"
    )

ax.set_xlabel("x (m)")
ax.set_ylabel("Temperature (K)")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("interface-temperature.png", dpi=300)
