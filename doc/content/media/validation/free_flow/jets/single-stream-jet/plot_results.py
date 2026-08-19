#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

os.chdir(SCRIPT_DIR)

D = 0.01

DATA = os.path.abspath(
    os.path.join(
        SCRIPT_DIR,
        "../../../../../../../validation/free_flow/jets/single-stream-jet/RTL10000",
    )
)
REF = os.path.join(DATA, "reference")
GOLD = os.path.join(DATA, "gold")

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"


def load_ref(name):
    return pd.read_csv(os.path.join(REF, name), header=None, names=["x", "y"])


def output_path(name):
    current = os.path.join(DATA, name)
    if os.path.exists(current):
        return current
    return os.path.join(GOLD, name)


def load_centerline_sim(path):
    df = pd.read_csv(path).sort_values("z").reset_index(drop=True)
    df["z_over_D"] = df["z"] / D
    df["Vin_over_Vc"] = df["vel_z"].iloc[0] / df["vel_z"]
    return df


def load_halfwidth_sim(path):
    df = pd.read_csv(path).sort_values("z").reset_index(drop=True)
    df["z_over_D"] = df["z"] / D
    df["rhalf_over_D"] = df["y_half"] / D
    return df


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

opgh_centerline = load_centerline_sim(
    output_path("jet_RTL10000_out_centerline_0002.csv")
)
opgh_halfwidth = load_halfwidth_sim(output_path("combined_halfwidths.csv"))

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
    label="OpenPronghorn (RealizableTwoLayer, Wolfstein)",
)

plt.plot(
    default_centerline["x"],
    default_centerline["y"],
    linestyle="-",
    linewidth=2.0,
    color="blue",
    label="Default IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    calibrated_centerline["x"],
    calibrated_centerline["y"],
    linestyle="-",
    linewidth=2.0,
    color="orange",
    label="Calibrated IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    lumley_centerline["x"],
    lumley_centerline["y"],
    linestyle="-",
    linewidth=2.0,
    color="green",
    label="Panchapakesan & Lumley 1993",
)

plt.plot(
    mi_centerline["x"],
    mi_centerline["y"],
    linestyle="None",
    marker="x",
    markersize=7,
    color="purple",
    label="Mi et al. 2013",
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
plt.close()

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
    label="Panchapakesan & Lumley 1993",
)

plt.plot(
    mi_halfwidth["x"],
    mi_halfwidth["y"],
    linestyle="None",
    marker="x",
    markersize=7,
    color="purple",
    label="Mi et al. 2013",
)

plt.plot(
    default_halfwidth["x"],
    default_halfwidth["y"],
    linestyle="-",
    linewidth=2.0,
    color="blue",
    label="Default IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    calibrated_halfwidth["x"],
    calibrated_halfwidth["y"],
    linestyle="-",
    linewidth=2.0,
    color="orange",
    label="Calibrated IP-RSM, Turutoglu et al. 2024",
)

plt.plot(
    opgh_halfwidth["z_over_D"],
    opgh_halfwidth["rhalf_over_D"],
    linestyle="-",
    linewidth=2.0,
    color="red",
    label="OpenPronghorn (RealizableTwoLayer, Wolfstein)",
)

plt.title(
    "Jet half-width growth\nSingle-Stream Jet in Still Air, Re = 10000",
    fontsize=13,
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
plt.close()
