import matplotlib.pyplot as plt
import pandas as pd
import os
from common import read_moose, analytic_approximate_solution_dht_erfcx, analytic_approximate_solution_bl

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Formal appearance
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.titlesize": 14,
    "axes.titleweight": 'bold'
})

# Read the files
print_commercial_cfd = False
moose_fluid_df = read_moose("../../../../../../../validation/free_flow/heat_transfer/flat_plate/gold/cht_rob-rob_out_y_vs_tf_0001.csv", "y")
moose_solid_df = read_moose("../../../../../../../validation/free_flow/heat_transfer/flat_plate/gold/cht_rob-rob_out_y_vs_ts_0001.csv", "y")
merged = pd.concat([moose_fluid_df, moose_solid_df])
merged = merged.sort_values("y")

merged['value_group'] = merged['T'].round(5)
merged = merged.groupby('value_group', as_index=False).mean()

if (print_commercial_cfd):
    reference_cfd_df = read_moose("../../../../../../../validation/free_flow/heat_transfer/flat_plate/gold/reference_cfd_horizontal.csv").sort_values("y")

# Get the analytic values for both approaches
st_bl = [analytic_approximate_solution_bl(0.1,y) for y in merged["y"]]
st_dht = [analytic_approximate_solution_dht_erfcx(0.1,y) for y in merged["y"]]

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(merged["y"], merged["T"], "b", label="MOOSE/Pronghorn")
ax.plot(merged["y"], st_bl, "r", label="Analytic BL (Luikov)")
ax.plot(merged["y"], st_dht, "g", label="Analytic DHT (Luikov)")

if (print_commercial_cfd):
    ax.plot(reference_cfd_df["y"], reference_cfd_df["T"], "--", label="Commercial CFD tool")

ax.set_xlabel("y (m)")
ax.set_ylabel("Temperature (K)")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig('vertical-temperature.png', dpi=300)
