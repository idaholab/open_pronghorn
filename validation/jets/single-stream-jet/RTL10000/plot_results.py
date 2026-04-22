#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

D = 0.01
REF = "reference"

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"


def load_ref(name):
    return pd.read_csv(os.path.join(REF, name), header=None, names=["x", "y"])


def load_centerline_sim(fname):
    df = pd.read_csv(fname).sort_values("z").reset_index(drop=True)
    df["z_over_D"] = df["z"] / D
    df["Vin_over_Vc"] = df["vel_z"].iloc[0] / df["vel_z"]
    return df


def load_halfwidth_sim(fname):
    df = pd.read_csv(fname).sort_values("z").reset_index(drop=True)
    df["z_over_D"] = df["z"] / D
    df["rhalf_over_D"] = df["y_half"] / D
    return df


# -------------------------
# Load reference data
# -------------------------
exp_centerline = load_ref("exp_centerline.csv")
default_centerline = load_ref("default_centerline.csv")
calibrated_centerline = load_ref("calibrated_centerline.csv")
balajee_centerline = load_ref("Balajee_centerline.csv")
lumley_centerline = load_ref("lumley_centerline.csv")
mi_centerline = load_ref("Mi_centerline.csv")

exp_halfwidth = load_ref("exp_halfwidth.csv")
lumley_halfwidth = load_ref("lumley_halfwidth.csv")
mi_halfwidth = load_ref("Mi_halfwidth.csv")
calibrated_halfwidth = load_ref("calibrated_IP_RSM_halfwidth.csv")
default_halfwidth = load_ref("default_IP_RSM_halfwidth.csv")

# -------------------------
# Load simulation data
# -------------------------
opgh_centerline = load_centerline_sim("jet_RTL10000_out_centerline_0002.csv")
opgh_halfwidth = load_halfwidth_sim("combined_halfwidths.csv")

# -------------------------
# Centerline plot
# -------------------------
plt.figure()

plt.plot(
    exp_centerline["x"],
    exp_centerline["y"],
    linestyle="None",
    marker="o",
    markersize=6,
    markerfacecolor="none",
    markeredgecolor="k",
    markeredgewidth=1.2,
    label="Experiment, Turutoglu et al. 2024",
)

plt.plot(
    opgh_centerline["z_over_D"],
    opgh_centerline["Vin_over_Vc"],
    linestyle="-",
    color="red",
    label="openPronghorn (RealizableTwoLayer, Wolfstein)",
)

plt.plot(
    default_centerline["x"],
    default_centerline["y"],
    linestyle="-",
    linewidth=2.0,
    color="blue",
    label="default IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    calibrated_centerline["x"],
    calibrated_centerline["y"],
    linestyle="-",
    linewidth=2.0,
    color="orange",
    label="calibrated IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    lumley_centerline["x"],
    lumley_centerline["y"],
    linestyle="-",
    linewidth=2.0,
    color="green",
    label="Panchapakesan & Lumley",
)

plt.plot(
    mi_centerline["x"],
    mi_centerline["y"],
    linestyle="None",
    marker="x",
    markersize=7,
    color="purple",
    label="Mi et al.",
)

plt.plot(
    balajee_centerline["x"],
    balajee_centerline["y"],
    linestyle="-",
    linewidth=2.0,
    color="pink",
    label="Balajee & Panchapakesan",
)

plt.title(
    "Center-line normalized axial velocity\nSingle-Stream Jet in Still Air, Re = 10000",
    fontsize=13,
)
plt.xlabel(r"$axial~location~[x/D]$", fontsize=14)
plt.ylabel(r"$V_{in}/V_x~[-]$", fontsize=14)
plt.legend(fontsize=8, loc="upper left")
plt.xlim(0, 25)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.savefig("jet_CL_normal_velocity.png", dpi=300, bbox_inches="tight")
plt.show()

# -------------------------
# Half-width plot
# -------------------------
plt.figure()

plt.plot(
    exp_halfwidth["x"],
    exp_halfwidth["y"],
    linestyle="None",
    marker="o",
    markersize=6,
    markerfacecolor="none",
    markeredgecolor="k",
    markeredgewidth=1.2,
    label="Experiment, Turutoglu et al. 2024",
)

plt.plot(
    lumley_halfwidth["x"],
    lumley_halfwidth["y"],
    linestyle="None",
    marker="+",
    markersize=7,
    color="green",
    label="Panchapakesan & Lumley correlation",
)

plt.plot(
    mi_halfwidth["x"],
    mi_halfwidth["y"],
    linestyle="None",
    marker="x",
    markersize=7,
    color="purple",
    label="Mi et al. correlation",
)

plt.plot(
    default_halfwidth["x"],
    default_halfwidth["y"],
    linestyle="-",
    linewidth=2.0,
    color="blue",
    label="default IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    calibrated_halfwidth["x"],
    calibrated_halfwidth["y"],
    linestyle="-",
    linewidth=2.0,
    color="orange",
    label="calibrated IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    opgh_halfwidth["z_over_D"],
    opgh_halfwidth["rhalf_over_D"],
    linestyle="-",
    linewidth=2.0,
    color="red",
    label="openPronghorn (RealizableTwoLayer, Wolfstein)",
)

plt.title(
    "Jet half-width growth\nSingle-Stream Jet in Still Air, Re = 10000", fontsize=13
)
plt.xlabel(r"$axial~location~[x/D]$", fontsize=14)
plt.ylabel(r"$r_{1/2}/D~[-]$", fontsize=14)
plt.legend(fontsize=8, loc="upper left")
plt.xlim(0, 25)
plt.ylim(0, 3)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.savefig("jet_halfwidth.png", dpi=300, bbox_inches="tight")
plt.show()
