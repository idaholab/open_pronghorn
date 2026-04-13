# %%
# This script plots the validation results based on the open-pronghorn simulation
# and experimental measurements. User needs to edit appropriately
##### LOAD MODULES ###############
import numpy as np
import matplotlib.pyplot as plt

##################################
##### Centerline velocity ########
exp_centerline = np.genfromtxt("exp_centerline.csv", skip_header=0, delimiter=",")
default_centerline = np.genfromtxt(
    "default_centerline.csv", skip_header=0, delimiter=","
)
calibrated_centerline = np.genfromtxt(
    "calibrated_centerline.csv", skip_header=0, delimiter=","
)
Balajee_centerline = np.genfromtxt(
    "Balajee_centerline.csv", skip_header=0, delimiter=","
)
lumley_centerline = np.genfromtxt("lumley_centerline.csv", skip_header=0, delimiter=",")
Mi_centerline = np.genfromtxt("Mi_centerline.csv", skip_header=0, delimiter=",")
opgh_RTL_centerline = np.genfromtxt(
    "../jet_RTL10000_out_centerline_0002.csv", skip_header=1, delimiter=","
)
#####
exp_halfwidth = np.genfromtxt("exp_halfwidth.csv", skip_header=0, delimiter=",")
opgh_RTL_halfwidth = np.genfromtxt(
    "../combined_halfwidths.csv", skip_header=1, delimiter=","
)
lumley_halfwidth = np.genfromtxt("lumley_halfwidth.csv", skip_header=0, delimiter=",")
Mi_halfwidth = np.genfromtxt("Mi_halfwidth.csv", skip_header=0, delimiter=",")
calibrated_IP_RSM_halfwidth = np.genfromtxt(
    "calibrated_IP_RSM_halfwidth.csv", skip_header=0, delimiter=","
)
default_IP_RSM_halfwidth = np.genfromtxt(
    "default_IP_RSM_halfwidth.csv", skip_header=0, delimiter=","
)
############### MAKE PRETTY ###################
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
###############################################

# plt.figure()

# --- Experimental data: Turutoglu, et al. 2024 ---
plt.plot(
    exp_centerline[:, 0],
    exp_centerline[:, 1],
    linestyle="None",
    marker="o",
    markersize=6,
    markerfacecolor="none",  # see-through
    markeredgecolor="k",  # black edge
    markeredgewidth=1.2,
    label="Experiment, Turutoglu et al. 2024",
)


# --- Simulation results RTL ---
plt.plot(
    opgh_RTL_centerline[:, -1] / 0.01,
    opgh_RTL_centerline[0, 1] / opgh_RTL_centerline[:, 1],
    linestyle="-",
    color="red",
    label="openPronghorn (RealizableTwoLayer, Wolfstein)",
    zorder=1,
)

# --- Default IP-RSM ---
plt.plot(
    default_centerline[:, 0],
    default_centerline[:, 1],
    linestyle="-",
    linewidth=2.0,
    color="blue",
    label="default IP-RSM, Turutoglu et al. 2024",
    zorder=1,
)

# --- Calibrated IP-RSM ---
plt.plot(
    calibrated_centerline[:, 0],
    calibrated_centerline[:, 1],
    linestyle="-",
    linewidth=2.0,
    color="orange",
    label="calibrated IP-RSM, Turutoglu et al. 2024",
    zorder=1,
)

# --- Lumley correlation ---
plt.plot(
    lumley_centerline[:, 0],
    lumley_centerline[:, 1],
    linestyle="-",
    linewidth=2.0,
    color="green",
    label="Panchapakesan & Lumley",
)

# --- Mi correlation ---
plt.plot(
    Mi_centerline[:, 0],
    Mi_centerline[:, 1],
    linestyle="None",
    marker="x",
    markersize=7,
    clip_on=False,
    color="purple",
    label="Mi et al.",
)

# --- Balajee ---
plt.plot(
    Balajee_centerline[:, 0],
    Balajee_centerline[:, 1],
    linestyle="-",
    linewidth=2.0,
    color="pink",
    label="Balajee & Panchapakesan",
)


# --- Titles and labels ---
plt.title(
    r"Center-line normalized axial velocity"
    "\n"
    "Single-Stream Jet in Still Air, Re = 10000",
    fontsize=13,
)
plt.xlabel(r"$axial~location~[x/D]$", fontsize=14)
plt.ylabel(r"$V_{in}/V_x~[-]$", fontsize=14)

# --- Legend and layout ---
plt.legend(fontsize=8, loc="upper left")
plt.xlim(0, 25)
# plt.ylim(0, 4.5)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)

# --- Save and show ---
plt.savefig("jet_CL_normal_velocity.png", dpi=300, bbox_inches="tight")
plt.show()

#########################################################################

plt.figure()

# --- Experimental data: Turutoglu, et al. 2024 ---
plt.plot(
    exp_halfwidth[:, 0],
    exp_halfwidth[:, 1],
    linestyle="None",
    marker="o",
    markersize=6,
    markerfacecolor="none",  # see-through
    markeredgecolor="k",  # black edge
    markeredgewidth=1.2,
    label="Experiment, Turutoglu et al. 2024",
)

# --- Lumley correlation ---
plt.plot(
    lumley_halfwidth[:, 0],
    lumley_halfwidth[:, 1],
    linestyle="None",
    marker="+",
    markersize=7,
    clip_on=False,
    color="green",
    label="Panchapakesan & Lumley correlation",
)

# --- Mi correlation ---
plt.plot(
    Mi_halfwidth[:, 0],
    Mi_halfwidth[:, 1],
    linestyle="None",
    marker="x",
    markersize=7,
    clip_on=False,
    color="purple",
    label="Mi et al. correlation",
)

# --- Default IP-RSM ---
plt.plot(
    default_IP_RSM_halfwidth[:, 0],
    default_IP_RSM_halfwidth[:, 1],
    linestyle="-",
    linewidth=2.0,
    color="blue",
    label="default IP-RSM, Turutoglu et al. 2024",
    zorder=1,
)

# --- Calibrated IP-RSM ---
plt.plot(
    calibrated_IP_RSM_halfwidth[:, 0],
    calibrated_IP_RSM_halfwidth[:, 1],
    linestyle="-",
    linewidth=2.0,
    color="orange",
    label="calibrated IP-RSM, Turutoglu et al. 2024",
    zorder=1,
)

# --- Pronghorn RTL ---
plt.plot(
    opgh_RTL_halfwidth[:, 2] / 0.01,
    opgh_RTL_halfwidth[:, 3] / 0.01,
    linestyle="-",
    linewidth=2.0,
    color="red",
    label="openPronghorn (RealizableTwoLayer, Wolfstein)",
    zorder=1,
)

# --- Titles and labels ---
plt.title(
    r"Jet half-width growth" "\n" "Single-Stream Jet in Still Air, Re = 10000",
    fontsize=13,
)
plt.xlabel(r"$axial~location~[x/D]$", fontsize=14)
plt.ylabel(r"$r_{1/2}/D~[-]$", fontsize=14)

# --- Legend and layout ---
plt.legend(fontsize=8, loc="upper left")
plt.xlim(0, 25)
plt.ylim(0, 3)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)

# --- Save and show ---
plt.savefig("jet_halfwidth.png", dpi=300, bbox_inches="tight")
plt.show()


# %%
