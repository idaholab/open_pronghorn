### ERCOFTAC Case 030: BFS with Straight Opposite Wall
### Plot Script
### Authors:
###  Dr. Mauricio Tanore
###  Hailey Tran-Kieu
###  Dr. Guillaume Giudicelli
### Date: October 24th, 2025

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# To make the relative paths work
os.chdir(os.path.dirname(os.path.realpath(__file__)))

plot_current_results = False
plot_reference_results = True
plot_errors = True

# Authorized relative increase in the error, shown in error bars
err_max = 0.01

############################ Model variants #############################
# Keys are internal, "base" is the filename stem in gold/, "label" is legend text
k_eps_variants = [
    {
        "key": "standard",
        "base": "k_epsilon_standard",
        "label": r"Std. k-$\varepsilon$",
    },
    {
        "key": "realizable",
        "base": "k_epsilon_realizable",
        "label": r"Realizable k-$\varepsilon$",
    },
    {
        "key": "realizable_yap",
        "base": "k_epsilon_realizable_yap",
        "label": r"Realizable k-$\varepsilon$ + YAP",
    },
    {
        "key": "realizable_two_layer",
        "base": "k_epsilon_realizable_two_layer",
        "label": r"Realizable k-$\varepsilon$ (two layer)",
    },
    {
        "key": "realizable_quadratic",
        "base": "k_epsilon_realizable_quadratic",
        "label": r"Realizable k-$\varepsilon$ (quadratic)",
    },
    {
        "key": "realizable_cubic",
        "base": "k_epsilon_realizable_cubic",
        "label": r"Realizable k-$\varepsilon$ (cubic)",
    },
]

############################ Plot settings #############################

# only add reference suffix when plotting the current results too
ref_suffix = ""
if plot_current_results:
    ref_suffix = " reference"

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

############################ Read .csv Data ############################

validation_folder = (
    "../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs/"
)

### Read simulation MOOSE .csv files (optional "current" run)
if plot_current_results:
    file_base = "bfs_input_csv"
    sim_in = pd.read_csv(
        validation_folder + file_base + "_inlet_channel_wall_sampler_0002.csv"
    )
    sim_out = pd.read_csv(
        validation_folder + file_base + "_outlet_channel_wall_sampler_0002.csv"
    )

### Read ALL reference MOOSE .csv files for each k-epsilon variant (wall samplers)
# Stored per model: x, pressure, mu_t_wall, distance, vel_x
wall_data = {}

for m in k_eps_variants:
    base = m["base"]
    inlet = pd.read_csv(
        validation_folder + f"gold/{base}_inlet_channel_wall_sampler_0002.csv"
    )
    outlet = pd.read_csv(
        validation_folder + f"gold/{base}_outlet_channel_wall_sampler_0002.csv"
    )

    x = np.concatenate([inlet["x"], outlet["x"]])
    pressure = np.concatenate([inlet["pressure"], outlet["pressure"]])
    mu_t = np.concatenate([inlet["mu_t_wall"], outlet["mu_t_wall"]])
    distance = np.concatenate([inlet["distance"], outlet["distance"]])
    vel_x = np.concatenate([inlet["vel_x"], outlet["vel_x"]])

    wall_data[m["key"]] = {
        "x": x,
        "pressure": pressure,
        "mu_t": mu_t,
        "distance": distance,
        "vel_x": vel_x,
    }

### Read ERCOFTAC benchmark .csv files
cp_exp = pd.read_csv(validation_folder + "reference_csv/cp.csv")
cf_exp = pd.read_csv(validation_folder + "reference_csv/cf.csv")

############################ Pressure Coefficient Analysis ############################

# Velocity at center of channel at x/H = -4
U_ref = 4.402663e01  # if your simulation is not exactly at 44.2 for free stream, update this number
# Height of the step
H = 0.0127
# Density of the fluid
rho = 1.18415
# Viscosity of the fluid
mu = 1.8551e-5
cp_factor = 0.5 * rho * U_ref**2  # 1/2 rho U_ref^2

# Uncertainty reported in http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case030
uncertainty_cp = 0.009

if plot_current_results:
    sim_x = np.concatenate([sim_in["x"], sim_out["x"]])
    sim_pressure = np.concatenate([sim_in["pressure"], sim_out["pressure"]])

# Create subplots
fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))

ax1.plot(cp_exp["x/h"], cp_exp["cp"], "k.", label="ERCOFTAC")
ax1.errorbar(
    cp_exp["x/h"],
    cp_exp["cp"],
    yerr=uncertainty_cp,
    ecolor="r",
    capsize=3,
    fmt="k.",
    barsabove=False,
    label="Reported 95% confidence",
)

# Plot each k-epsilon variant
if plot_reference_results:
    for m in k_eps_variants:
        d = wall_data[m["key"]]
        x = d["x"]
        pressure = d["pressure"]
        cp_model = pressure / cp_factor + 0.125
        ax1.plot(
            x / H,
            cp_model,
            label=f"OpenPronghorn {m['label']}{ref_suffix}",
        )

if plot_current_results:
    ax1.plot(
        sim_x / H,
        sim_pressure / cp_factor + 0.125,
        "b",
        label="OpenPronghorn k-epsilon current",
    )

ax1.set_xlabel(r"$\mathrm{x/h}$", fontsize=14)
ax1.set_ylabel(r"$\mathrm{c_p}$", fontsize=14)
ax1.set_xlim([-5, 22])
ax1.grid(True)

# Legend outside (to the right)
ax1.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

plt.title("Pressure Coefficient")
fig.tight_layout()
fig.savefig("plots_cp_main.png", dpi=300, bbox_inches="tight")

############################ Wall-Skin Friction Coefficient Analysis ############################

cf_factor = 0.5 * rho * U_ref**2

if plot_current_results:
    sim_x = np.concatenate([sim_in["x"], sim_out["x"]])
    sim_mu_t = np.concatenate([sim_in["mu_t_wall"], sim_out["mu_t_wall"]])
    sim_distance = np.concatenate([sim_in["distance"], sim_out["distance"]])
    sim_vel_x = np.concatenate([sim_in["vel_x"], sim_out["vel_x"]])

# Uncertainty reported in http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case030
uncertainty_cf_inlet_percent = 0.08
uncertainty_cf_detached_percent = 0.15
uncertainty_cf = np.abs(cf_exp["cf"]) * np.array(
    [
        uncertainty_cf_inlet_percent * (x_over_H < 0 or x_over_H > 6.26)
        + uncertainty_cf_detached_percent * (x_over_H >= 0 and x_over_H <= 6.26)
        for x_over_H in cf_exp["x/h"]
    ]
)

# Create subplots
fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))

ax1.plot(cf_exp["x/h"], cf_exp["cf"], "k.", label="Exp.")
ax1.errorbar(
    cf_exp["x/h"],
    cf_exp["cf"],
    yerr=uncertainty_cf,
    ecolor="r",
    capsize=3,
    fmt="k.",
    barsabove=False,
    label="Reported 95% confidence",
)

# Plot each k-epsilon variant
if plot_reference_results:
    for m in k_eps_variants:
        d = wall_data[m["key"]]
        x = d["x"]
        mu_t = d["mu_t"]
        distance = d["distance"]
        vel_x = d["vel_x"]
        cf_model = ((mu_t + mu) * vel_x / distance) / cf_factor

        ax1.plot(
            x / H,
            cf_model,
            linestyle="-",
            marker=None,
            label=f"OpenPronghorn {m['label']}{ref_suffix}",
        )

if plot_current_results:
    ax1.plot(
        sim_x / H,
        ((sim_mu_t + mu) * sim_vel_x / sim_distance) / cf_factor,
        "b",
        label="OpenPronghorn k-epsilon current",
    )

ax1.set_xlabel(r"$\mathrm{x/h}$", fontsize=14)
ax1.set_ylabel(r"$\mathrm{c_f}$", fontsize=14)
ax1.set_xlim([-5, 22])
ax1.grid(True)

# Legend outside
ax1.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)

plt.title("Skin Friction Coefficient")
fig.tight_layout()
fig.savefig("plots_cf_main.png", dpi=300, bbox_inches="tight")

############################ Vertical X-Velocity Profiles Analysis ############################

# Read the y_grid values
y_grid = pd.read_csv(validation_folder + "reference_csv/u+1.csv")["y/h"].to_numpy()

# MOOSE Reference: Function to interpolate vel_x to y_grid
def interpolate_vel_x(y_values, vel_x_values, y_grid):
    return np.interp(y_grid, y_values, vel_x_values)


# For each model, read vel_x profiles at x/H = 1,4,6,10 and interpolate
x_locations = [1, 4, 6, 10]
vel_profiles = {m["key"]: {} for m in k_eps_variants}

for m in k_eps_variants:
    base = m["base"]
    for xh in x_locations:
        file = f"gold/{base}_vel_x_xoH_{xh}_sampler_0002.csv"
        df_linear = pd.read_csv(validation_folder + file)
        y_linear = (df_linear["y"].to_numpy() + H) / H  # Normalize y-values to y/H
        vel_x_linear = df_linear["vel_x"].to_numpy() / U_ref  # Normalize velocities
        vel_profiles[m["key"]][xh] = interpolate_vel_x(
            y_linear, vel_x_linear, y_grid
        )

# Read and process each linear CSV for the current run
sim_vel_profiles = {}
if plot_current_results:
    sim_vel_profiles = {}
    sim_linear_files = {
        xh: f"{file_base}_vel_x_xoH_{xh}_sampler_0002.csv" for xh in x_locations
    }
    for xh, file in sim_linear_files.items():
        sim_df_linear = pd.read_csv(validation_folder + file)
        sim_y_linear = (sim_df_linear["y"].to_numpy() + H) / H
        sim_vel_x_linear = sim_df_linear["vel_x"].to_numpy() / U_ref
        sim_vel_profiles[xh] = interpolate_vel_x(
            sim_y_linear, sim_vel_x_linear, y_grid
        )

# Read u profiles
u_exp = pd.read_csv(validation_folder + "reference_csv/u_profiles_exp.csv")

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    2, 2, figsize=(10, 10)
)
fig.suptitle("Vertical X-Velocity Profiles", fontsize="14")

# Map x/H to exp column
exp_cols = {1: "x+1", 4: "x+4", 6: "x+6", 10: "x+10"}
axes = {1: ax1, 4: ax2, 6: ax3, 10: ax4}

for xh in x_locations:
    ax = axes[xh]
    col = exp_cols[xh]

    ax.plot(u_exp[col], u_exp["y/h"], "k.", label="ERCOFTAC")

    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][xh]
            ax.plot(
                u_model,
                y_grid,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )

    if plot_current_results and xh in sim_vel_profiles:
        ax.plot(
            sim_vel_profiles[xh],
            y_grid,
            "b",
            label="OpenPronghorn k-epsilon current",
        )

    ax.set_title(f"Vertical x-Velocity Profile for x/h = {xh}", fontsize=14)
    ax.set_xlabel(r"$\mathrm{U/U_{ref}}$", fontsize=14)
    ax.set_ylabel(r"$\mathrm{y/h}$", fontsize=14)
    ax.set_ylim([-0.1, 4.5])
    ax.grid(True)

# Build a combined legend outside (to the right)
handles, labels = ax1.get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    borderaxespad=0.0,
    fontsize=8,
)

fig.tight_layout()
fig.savefig("plots_u_profiles_main.png", dpi=300, bbox_inches="tight")

############################ Generating Error Bar Plots ############################

############################ Pressure Coefficient and Skin Friction Coefficient ############################

ercoftac_cp = cp_exp["cp"].to_numpy()
ercoftac_cf = cf_exp["cf"].to_numpy()
ercoftac_x_cp = cp_exp["x/h"].to_numpy()
ercoftac_x_cf = cf_exp["x/h"].to_numpy()

# Interpolate each model onto ERCOFTAC x locations (cp and cf)
cp_errors = {}  # key -> array (model - exp)
cf_errors = {}

if plot_reference_results:
    for m in k_eps_variants:
        d = wall_data[m["key"]]
        moose_x = d["x"] / H

        moose_cp = d["pressure"] / cp_factor + 0.125
        moose_cf = ((d["mu_t"] + mu) * d["vel_x"] / d["distance"]) / cf_factor

        moose_cp_interp = np.interp(ercoftac_x_cp, moose_x, moose_cp)
        moose_cf_interp = np.interp(ercoftac_x_cf, moose_x, moose_cf)

        cp_errors[m["key"]] = moose_cp_interp - ercoftac_cp
        cf_errors[m["key"]] = moose_cf_interp - ercoftac_cf

# Current run error (optional)
if plot_current_results:
    sim_moose_x = sim_x / H
    sim_moose_cp = sim_pressure / cp_factor + 0.125
    sim_moose_cf = ((sim_mu_t + mu) * sim_vel_x / sim_distance) / cf_factor

    sim_moose_cp_interp = np.interp(ercoftac_x_cp, sim_moose_x, sim_moose_cp)
    sim_moose_cf_interp = np.interp(ercoftac_x_cf, sim_moose_x, sim_moose_cf)

fig, (ax5, ax6) = plt.subplots(1, 2, figsize=(10, 5))

if not plot_errors:
    # Use first variant as reference for error-band construction
    ref_key = k_eps_variants[0]["key"]
    moose_cp_interp_ref = cp_errors[ref_key] + ercoftac_cp
    moose_cf_interp_ref = cf_errors[ref_key] + ercoftac_cf

    error_magnitude_cp = np.abs((1 + err_max) * (ercoftac_cp - moose_cp_interp_ref))
    error_magnitude_cf = np.abs((1 + err_max) * (ercoftac_cf - moose_cf_interp_ref))

    ax5.errorbar(
        ercoftac_x_cp,
        ercoftac_cp,
        yerr=error_magnitude_cp,
        ecolor="r",
        capsize=3,
        fmt="k-",
        barsabove=False,
    )
    ax5.plot(ercoftac_x_cp, ercoftac_cp, "k.", markersize=5, label="ERCOFTAC")
    if plot_reference_results:
        for m in k_eps_variants:
            ax5.plot(
                ercoftac_x_cp,
                cp_errors[m["key"]] + ercoftac_cp,
                linestyle="--",
                marker="o",
                markersize=3,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results:
        ax5.plot(
            ercoftac_x_cp,
            sim_moose_cp_interp,
            color="b",
            linestyle="--",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
    ax5.set_ylabel(r"$\mathrm{c_p}$", fontsize=14)
    ax5.set_title("Pressure Coefficient with Error Bars")

    ax6.errorbar(
        ercoftac_x_cf,
        ercoftac_cf,
        yerr=error_magnitude_cf,
        ecolor="r",
        capsize=3,
        fmt="k-",
        barsabove=False,
    )
    ax6.plot(ercoftac_x_cf, ercoftac_cf, "k.", markersize=5, label="ERCOFTAC")
    if plot_reference_results:
        for m in k_eps_variants:
            ax6.plot(
                ercoftac_x_cf,
                cf_errors[m["key"]] + ercoftac_cf,
                linestyle="--",
                marker="o",
                markersize=3,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results:
        ax6.plot(
            ercoftac_x_cf,
            sim_moose_cf_interp,
            color="b",
            linestyle="--",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
    ax6.set_ylabel(r"$\mathrm{c_f}$", fontsize=14)
    ax6.set_title("Wall-Skin Friction Coefficient with Error Bars")

else:
    # Pure error plots: model - experimental
    if plot_reference_results:
        for m in k_eps_variants:
            ax5.plot(
                ercoftac_x_cp,
                cp_errors[m["key"]],
                linestyle="--",
                marker="o",
                markersize=3,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
            ax6.plot(
                ercoftac_x_cf,
                cf_errors[m["key"]],
                linestyle="--",
                marker="o",
                markersize=3,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )

    if plot_current_results:
        ax5.plot(
            ercoftac_x_cp,
            sim_moose_cp_interp - ercoftac_cp,
            color="b",
            linestyle="--",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
        ax6.plot(
            ercoftac_x_cf,
            sim_moose_cf_interp - ercoftac_cf,
            color="b",
            linestyle="--",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )

    ax5.set_ylabel(r"$\mathrm{c_p - c_{p,exp}}$", fontsize=14)
    ax5.set_title("Error in Pressure Coefficient")
    ax6.set_ylabel(r"$\mathrm{c_f - c_{f,exp}}$", fontsize=14)
    ax6.set_title("Error in Wall-Skin Friction Coefficient")

ax5.set_xlabel(r"$\mathrm{x/h}$", fontsize=14)
ax6.set_xlabel(r"$\mathrm{x/h}$", fontsize=14)
ax5.set_xlim([-5, 22])
ax6.set_xlim([-5, 22])
ax5.grid(True)
ax6.grid(True)

# Legend outside (shared, to the right)
handles, labels = ax5.get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    borderaxespad=0.0,
    fontsize=8,
)

fig.tight_layout()
fig.savefig("plots_cp_cf_error_main.png", dpi=300, bbox_inches="tight")

############################ Velocity Profile Error Plots ############################

# Extract necessary experimental data
u_exp_1 = u_exp["x+1"].to_numpy()
u_exp_4 = u_exp["x+4"].to_numpy()
u_exp_6 = u_exp["x+6"].to_numpy()
u_exp_10 = u_exp["x+10"].to_numpy()

# Clean NaNs (assuming at most trailing NaNs; y_grid still used for vertical axis)
u_exp_1_cleaned = u_exp_1[~np.isnan(u_exp_1)]
u_exp_4_cleaned = u_exp_4[~np.isnan(u_exp_4)]
u_exp_6_cleaned = u_exp_6[~np.isnan(u_exp_6)]
u_exp_10_cleaned = u_exp_10[~np.isnan(u_exp_10)]

# Choose a reference variant for non-error error bars
ref_key = k_eps_variants[0]["key"]
u_ref_1 = vel_profiles[ref_key][1]
u_ref_4 = vel_profiles[ref_key][4]
u_ref_6 = vel_profiles[ref_key][6]
u_ref_10 = vel_profiles[ref_key][10]

# Align lengths crudely (assuming same vertical grid except for NaNs)
error_magnitude_1 = np.abs(
    (1 + err_max) * (u_exp_1_cleaned - u_ref_1[: len(u_exp_1_cleaned)])
)
error_magnitude_4 = np.abs(
    (1 + err_max) * (u_exp_4_cleaned - u_ref_4[: len(u_exp_4_cleaned)])
)
error_magnitude_6 = np.abs(
    (1 + err_max) * (u_exp_6_cleaned - u_ref_6[: len(u_exp_6_cleaned)])
)
error_magnitude_10 = np.abs(
    (1 + err_max) * (u_exp_10_cleaned - u_ref_10[: len(u_exp_10_cleaned)])
)

if plot_errors:
    x_label = r"$\mathrm{(U - U_{exp})/U_{ref}}$"
else:
    x_label = r"$\mathrm{U/U_{ref}}$"

# x/H = 1 and x/H = 4
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
if plot_errors:
    fig.suptitle("Error in Vertical X-Velocity Profiles", fontsize="14")
else:
    fig.suptitle("Vertical X-Velocity Profiles with Error Bars", fontsize="14")

# x/H = 1
if not plot_errors:
    ax1.errorbar(
        u_exp_1_cleaned,
        y_grid[: len(u_exp_1_cleaned)],
        xerr=error_magnitude_1,
        ecolor="r",
        capsize=3,
        fmt="k.",
        barsabove=False,
    )
    ax1.plot(
        u_exp_1_cleaned,
        y_grid[: len(u_exp_1_cleaned)],
        marker="o",
        markersize=4,
        label="ERCOFTAC",
    )
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][1][: len(u_exp_1_cleaned)]
            ax1.plot(
                u_model,
                y_grid[: len(u_exp_1_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 1 in sim_vel_profiles:
        sim_u = sim_vel_profiles[1][: len(u_exp_1_cleaned)]
        ax1.plot(
            sim_u,
            y_grid[: len(u_exp_1_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
else:
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][1][: len(u_exp_1_cleaned)]
            ax1.plot(
                u_model - u_exp_1_cleaned,
                y_grid[: len(u_exp_1_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 1 in sim_vel_profiles:
        sim_u = sim_vel_profiles[1][: len(u_exp_1_cleaned)]
        ax1.plot(
            sim_u - u_exp_1_cleaned,
            y_grid[: len(u_exp_1_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
ax1.set_xlabel(x_label, fontsize=14)
ax1.set_ylabel(r"$\mathrm{y/h}$", fontsize=14)
ax1.grid(True)
ax1.set_title("x/H = 1")
ax1.set_ylim([-0.1, 4.5])

# x/H = 4
if not plot_errors:
    ax2.errorbar(
        u_exp_4_cleaned,
        y_grid[: len(u_exp_4_cleaned)],
        xerr=error_magnitude_4,
        ecolor="r",
        capsize=3,
        fmt="k.",
        barsabove=False,
    )
    ax2.plot(
        u_exp_4_cleaned,
        y_grid[: len(u_exp_4_cleaned)],
        marker="o",
        markersize=4,
        label="ERCOFTAC",
    )
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][4][: len(u_exp_4_cleaned)]
            ax2.plot(
                u_model,
                y_grid[: len(u_exp_4_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 4 in sim_vel_profiles:
        sim_u = sim_vel_profiles[4][: len(u_exp_4_cleaned)]
        ax2.plot(
            sim_u,
            y_grid[: len(u_exp_4_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
else:
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][4][: len(u_exp_4_cleaned)]
            ax2.plot(
                u_model - u_exp_4_cleaned,
                y_grid[: len(u_exp_4_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 4 in sim_vel_profiles:
        sim_u = sim_vel_profiles[4][: len(u_exp_4_cleaned)]
        ax2.plot(
            sim_u - u_exp_4_cleaned,
            y_grid[: len(u_exp_4_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
ax2.set_xlabel(x_label, fontsize=14)
ax2.set_ylabel(r"$\mathrm{y/h}$", fontsize=14)
ax2.grid(True)
ax2.set_title("x/H = 4")
ax2.set_ylim([-0.1, 4.5])

# Legend outside (shared)
handles, labels = ax1.get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    borderaxespad=0.0,
    fontsize=8,
)

fig.tight_layout()
fig.savefig("plots_u_profiles_error1_main.png", dpi=300, bbox_inches="tight")

# x/H = 6 and x/H = 10
fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 5))
if plot_errors:
    fig.suptitle("Error in Vertical X-Velocity Profiles", fontsize="14")
else:
    fig.suptitle("Vertical X-Velocity Profiles with Error Bars", fontsize="14")

# x/H = 6
if not plot_errors:
    ax3.errorbar(
        u_exp_6_cleaned,
        y_grid[: len(u_exp_6_cleaned)],
        xerr=error_magnitude_6,
        ecolor="r",
        capsize=3,
        fmt="k.",
        barsabove=False,
    )
    ax3.plot(
        u_exp_6_cleaned,
        y_grid[: len(u_exp_6_cleaned)],
        marker="o",
        markersize=4,
        label="ERCOFTAC",
    )
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][6][: len(u_exp_6_cleaned)]
            ax3.plot(
                u_model,
                y_grid[: len(u_exp_6_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 6 in sim_vel_profiles:
        sim_u = sim_vel_profiles[6][: len(u_exp_6_cleaned)]
        ax3.plot(
            sim_u,
            y_grid[: len(u_exp_6_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
else:
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][6][: len(u_exp_6_cleaned)]
            ax3.plot(
                u_model - u_exp_6_cleaned,
                y_grid[: len(u_exp_6_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 6 in sim_vel_profiles:
        sim_u = sim_vel_profiles[6][: len(u_exp_6_cleaned)]
        ax3.plot(
            sim_u - u_exp_6_cleaned,
            y_grid[: len(u_exp_6_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
ax3.set_xlabel(x_label, fontsize=14)
ax3.set_ylabel(r"$\mathrm{y/h}$", fontsize=14)
ax3.grid(True)
ax3.set_title("x/H = 6")
ax3.set_ylim([-0.1, 4.5])

# x/H = 10
if not plot_errors:
    ax4.errorbar(
        u_exp_10_cleaned,
        y_grid[: len(u_exp_10_cleaned)],
        xerr=error_magnitude_10,
        ecolor="r",
        capsize=3,
        fmt="k.",
        barsabove=False,
    )
    ax4.plot(
        u_exp_10_cleaned,
        y_grid[: len(u_exp_10_cleaned)],
        marker="o",
        markersize=4,
        label="ERCOFTAC",
    )
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][10][: len(u_exp_10_cleaned)]
            ax4.plot(
                u_model,
                y_grid[: len(u_exp_10_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 10 in sim_vel_profiles:
        sim_u = sim_vel_profiles[10][: len(u_exp_10_cleaned)]
        ax4.plot(
            sim_u,
            y_grid[: len(u_exp_10_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
else:
    if plot_reference_results:
        for m in k_eps_variants:
            u_model = vel_profiles[m["key"]][10][: len(u_exp_10_cleaned)]
            ax4.plot(
                u_model - u_exp_10_cleaned,
                y_grid[: len(u_exp_10_cleaned)],
                linestyle="--",
                marker="o",
                markersize=4,
                label=f"OpenPronghorn {m['label']}{ref_suffix}",
            )
    if plot_current_results and 10 in sim_vel_profiles:
        sim_u = sim_vel_profiles[10][: len(u_exp_10_cleaned)]
        ax4.plot(
            sim_u - u_exp_10_cleaned,
            y_grid[: len(u_exp_10_cleaned)],
            linestyle="--",
            color="b",
            marker="o",
            markersize=4,
            label="OpenPronghorn k-epsilon current",
        )
ax4.set_xlabel(x_label, fontsize=14)
ax4.set_ylabel(r"$\mathrm{y/h}$", fontsize=14)
ax4.grid(True)
ax4.set_title("x/H = 10")
ax4.set_ylim([-0.1, 4.5])

# Legend outside (shared)
handles, labels = ax3.get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    borderaxespad=0.0,
    fontsize=8,
)

fig.tight_layout()
fig.savefig("plots_u_profiles_error2_main.png", dpi=300, bbox_inches="tight")
