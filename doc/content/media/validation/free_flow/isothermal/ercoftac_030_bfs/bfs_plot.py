### ERCOFTAC Case 030: BFS with Straight Opposite Wall
### Plot Script courtesy of Dr. Mauricio Tanoret
### Last Edited by Hailey Tran-Kieu
### Date: August 12, 2025

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# To make the relative paths work
os.chdir(os.path.dirname(os.path.realpath(__file__)))

plot_current_results = False
plot_reference_results = True

# Authorized relative increase in the error
err_max = 0.01

############################ Plot settings #############################

# only add reference suffix when plotting the current results too
ref_suffix = ''
if plot_current_results:
    ref_suffix = ' reference'

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

############################ Read .csv Data ############################

validation_folder = '../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs/'

### Read simulation MOOSE .csv files
# Modify based on current MOOSE results
if plot_current_results:
    sim_in = pd.read_csv(validation_folder + 'bfs_input_csv_inlet_sampler_0002.csv')
    sim_out = pd.read_csv(validation_folder + 'bfs_input_csv_outlet_sampler_0002.csv')

### Read reference MOOSE .csv files
# Modify if changing reference parameters
reference_in = pd.read_csv(validation_folder + 'reference_csv/reference_bfs_input_csv_inlet_sampler_0002.csv')
reference_out = pd.read_csv(validation_folder + 'reference_csv/reference_bfs_input_csv_outlet_sampler_0002.csv')

### Read ERCOFTAC benchmark .csv files
cp_exp = pd.read_csv(validation_folder + 'reference_csv/cp.csv')
cf_exp = pd.read_csv(validation_folder + 'reference_csv/cf.csv')

############################ Pressure Coefficient Analysis ############################

H = 0.0127
cp_factor = 0.5*1.18415*48.18**2 # 1/2 rho U_ref^2

# Concatenate inlet and outlet data for MOOSE .csv
x = np.concatenate([reference_in['x'], reference_out['x']])
pressure = np.concatenate([reference_in['pressure'], reference_out['pressure']])

if plot_current_results:
    sim_x = np.concatenate([sim_in['x'], sim_out['x']])
    sim_pressure = np.concatenate([sim_in['pressure'], sim_out['pressure']])


# Create subplots
fig, ax1 = plt.subplots(1, 1, figsize=(4, 4))

ax1.plot(cp_exp['x/h'], cp_exp['cp'], 'k.', label='ERCOFTAC')
if plot_reference_results:
    ax1.plot(x/H, pressure/cp_factor + 0.125, 'g', label='MOOSE-FV k-epsilon' + ref_suffix)
if plot_current_results:
    ax1.plot(sim_x/H, sim_pressure/cp_factor + 0.125, 'b', label='MOOSE-FV k-epsilon')
ax1.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax1.set_xlim([-5,22])
ax1.legend()
ax1.grid(True)

plt.title('Pressure Coefficient')
plt.tight_layout()
plt.savefig('plots_cp_main.png', dpi=300)

############################ Wall-Skin Friction Coefficient Analysis ############################

H = 0.0127
cf_factor = 0.5*1.18415*47.**2


# Concatenate inlet and outlet data for MOOSE .csv
x = np.concatenate([reference_in['x'], reference_out['x']])
mu_t = np.concatenate([reference_in['mu_t'], reference_out['mu_t']]) / 6.0
distance = np.concatenate([reference_in['distance'], reference_out['distance']])
vel_x = np.concatenate([reference_in['vel_x'], reference_out['vel_x']])

if plot_current_results:
    sim_x = np.concatenate([sim_in['x'], sim_out['x']])
    sim_mu_t = np.concatenate([sim_in['mu_t'], sim_out['mu_t']]) / 6.0
    sim_distance = np.concatenate([sim_in['distance'], sim_out['distance']])
    sim_vel_x = np.concatenate([sim_in['vel_x'], sim_out['vel_x']])


# Create subplots
fig, ax1 = plt.subplots(1, 1, figsize=(4, 4))

ax1.plot(cf_exp['x/h'], cf_exp['cf'], 'k.', label='Exp')
if plot_reference_results:
    ax1.plot(x/H, (mu_t*vel_x/distance)/cf_factor, 'g', label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax1.plot(sim_x/H, (sim_mu_t*sim_vel_x/sim_distance)/cf_factor, 'b', label='MOOSE-FV K-epsilon current')

ax1.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{c_f}$', fontsize=14)
ax1.set_xlim([-5,22])
ax1.legend()
ax1.grid(True)

# Adjust layout
plt.title('Skin Friction Coefficient')
plt.tight_layout()
plt.savefig('plots_cf_main.png', dpi=300)





############################ Vertical X-Velocity Profiles Analysis ############################

# Read the y_grid values
y_grid = pd.read_csv(validation_folder + 'reference_csv/u+1.csv')['y/h'].to_numpy()

# MOOSE Reference: Function to interpolate vel_x to y_grid
def interpolate_vel_x(y_values, vel_x_values, y_grid):
    return np.interp(y_grid, y_values, vel_x_values)

# Read and process each linear CSV
linear_files = [
    'reference_csv/reference_bfs_input_csv_vel_x_xoH_1_sampler_0002.csv',
    'reference_csv/reference_bfs_input_csv_vel_x_xoH_4_sampler_0002.csv',
    'reference_csv/reference_bfs_input_csv_vel_x_xoH_6_sampler_0002.csv',
    'reference_csv/reference_bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'
]

interpolated_results = {}

for file in linear_files:
    df_linear = pd.read_csv(validation_folder + file)
    y_linear = (df_linear['y'].to_numpy() + 0.0127) / 0.0127 # Normalize y-values to y/H
    vel_x_linear = df_linear['vel_x'].to_numpy() / 48.2 # Normalize velocities
    interpolated_vel_x = interpolate_vel_x(y_linear, vel_x_linear, y_grid) # Interpolate MOOSE data to y/H
    interpolated_results[file] = interpolated_vel_x


# Read and process each linear CSV
sim_linear_files = [
    'bfs_input_csv_vel_x_xoH_1_sampler_0002.csv',
    'bfs_input_csv_vel_x_xoH_4_sampler_0002.csv',
    'bfs_input_csv_vel_x_xoH_6_sampler_0002.csv',
    'bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'
]

sim_interpolated_results = {}

if plot_current_results:
    for file in sim_linear_files:
        sim_df_linear = pd.read_csv(validation_folder + file)
        sim_y_linear = (sim_df_linear['y'].to_numpy() + 0.0127) / 0.0127 # Normalize y-values to y/H
        sim_vel_x_linear = sim_df_linear['vel_x'].to_numpy() / 48.2 # Normalize velocities
        sim_interpolated_vel_x = interpolate_vel_x(sim_y_linear, sim_vel_x_linear, y_grid) # Interpolate MOOSE data to y/H
        sim_interpolated_results[file] = sim_interpolated_vel_x


# Read u profiles
u_exp = pd.read_csv(validation_folder + 'reference_csv/u_profiles_exp.csv')


# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(4, 4), constrained_layout=True)
fig.suptitle("Vertical X-Velocity Profiles", fontsize="14")

ax1.plot(u_exp['x+1'], u_exp['y/h'], 'k.', label='ERCOFTAC')
if plot_reference_results:
    ax1.plot(interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_1_sampler_0002.csv'], y_grid, 'g', label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax1.plot(sim_interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv'], y_grid, 'b', label='MOOSE-FV K-epsilon current')
ax1.set_title('Vertical x-Velocity Profile for x/h = 1', fontsize=14)
ax1.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax1.set_ylim([-0.1, 4.5])
ax1.legend(fontsize="10")
ax1.grid(True)

ax2.plot(u_exp['x+4'], u_exp['y/h'], 'k.', label='ERCOFTAC')
if plot_reference_results:
    ax2.plot(interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_4_sampler_0002.csv'], y_grid, 'g', label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax2.plot(sim_interpolated_results['bfs_input_csv_vel_x_xoH_4_sampler_0002.csv'], y_grid, 'b', label='MOOSE-FV K-epsilon current')
ax2.set_title('Vertical x-Velocity Profile for x/h = 4', fontsize=14)
ax2.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax2.set_ylim([-0.1, 4.5])
ax2.legend(fontsize="10")
ax2.grid(True)

ax3.plot(u_exp['x+6'], u_exp['y/h'], 'k.', label='ERCOFTAC')
if plot_reference_results:
    ax3.plot(interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_6_sampler_0002.csv'], y_grid, 'g', label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax3.plot(sim_interpolated_results['bfs_input_csv_vel_x_xoH_6_sampler_0002.csv'], y_grid, 'b', label='MOOSE-FV K-epsilon current')
ax3.set_title('Vertical x-Velocity Profile for x/h = 6', fontsize=14)
ax3.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax3.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax3.set_ylim([-0.1, 4.5])
ax3.legend(fontsize="10")
ax3.grid(True)

ax4.plot(u_exp['x+10'], u_exp['y/h'], 'k.', label='ERCOFTAC')
if plot_reference_results:
    ax4.plot(interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'], y_grid, 'g', label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax4.plot(sim_interpolated_results['bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'], y_grid, 'b', label='MOOSE-FV K-epsilon current')
ax4.set_title('Vertical x-Velocity Profile for x/h = 10', fontsize=14)
ax4.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax4.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax4.set_ylim([-0.1, 4.5])
ax4.legend(fontsize="10")
ax4.grid(True)


# Adjust layout
plt.tight_layout()
plt.savefig('plots_u_profiles_main.png', dpi=300)


############################ Generating Error Bar Plots ############################

############################ Pressure Coefficient and Skin Friction Coefficient ############################
ercoftac_cp = cp_exp['cp']
moose_cp = (pressure / cp_factor) + 0.125
sim_moose_cp = (pressure / cp_factor) + 0.125

ercoftac_cf = cf_exp['cf']
moose_cf = (mu_t*vel_x/distance)/cf_factor
sim_moose_cf = (mu_t*vel_x/distance)/cf_factor

moose_x = x/H
sim_moose_x = x/H
ercoftac_x_cp = cp_exp['x/h']
ercoftac_x_cf = cf_exp['x/h']

### Interpolate MOOSE onto ERCOFTAC
sim_moose_cp_interp = np.interp(cp_exp['x/h'], sim_moose_x, sim_moose_cp)
moose_cp_interp = np.interp(cp_exp['x/h'], moose_x, moose_cp)
error_magnitude_cp = (abs( (1 + err_max) * (ercoftac_cp - moose_cp_interp) ) )

sim_moose_cf_interp = np.interp(cf_exp['x/h'], sim_moose_x, sim_moose_cf)
moose_cf_interp = np.interp(cf_exp['x/h'], moose_x, moose_cf)
error_magnitude_cf = (abs( (1 + err_max) * (ercoftac_cf - moose_cf_interp) ) )

# Error Ranges
min_error_cp = ercoftac_cp - error_magnitude_cp
max_error_cp = ercoftac_cp + error_magnitude_cp

min_error_cf = ercoftac_cf - error_magnitude_cf
max_error_cf = ercoftac_cf + error_magnitude_cf

# Create subplots
fig, (ax5, ax6) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax5.errorbar(ercoftac_x_cp, ercoftac_cp, yerr=error_magnitude_cp, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax5.plot(ercoftac_x_cp, ercoftac_cp, 'k.', markersize=5, label='ERCOFTAC')
if plot_reference_results:
    ax5.plot(ercoftac_x_cp, moose_cp_interp, color='g', linestyle='--', marker='o', markersize=4, label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax5.plot(ercoftac_x_cp, sim_moose_cp_interp, color='b', linestyle='--', marker='o', markersize=4, label='MOOSE-FV K-epsilon current')
ax5.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax5.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax5.legend()
ax5.grid(True)
ax5.set_title("Pressure Coefficient with Error Bars")
ax5.set_xlim([-5,22])

ax6.errorbar(ercoftac_x_cf, ercoftac_cf, yerr=error_magnitude_cf, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax6.plot(ercoftac_x_cf, ercoftac_cf, 'k.', markersize=5, label='ERCOFTAC')
if plot_reference_results:
    ax6.plot(ercoftac_x_cf, moose_cf_interp, color='g', linestyle='--', marker='o', markersize=4, label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax6.plot(ercoftac_x_cf, sim_moose_cf_interp, color='b', linestyle='--', marker='o', markersize=4, label='MOOSE-FV K-epsilon current')
ax6.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax6.set_ylabel(r'$\mathrm{c_f}$', fontsize=14)
ax6.legend()
ax6.grid(True)
ax6.set_title("Wall-Skin Friction Coefficient with Error Bars")

plt.savefig('plots_cp_cf_error_main.png', dpi=300)


# Extract necessary data
u_exp_1 = u_exp['x+1']
u_exp_1_cleaned = u_exp_1[~np.isnan(u_exp_1)]
u_exp_4 = u_exp['x+4']
u_exp_4_cleaned = u_exp_4[~np.isnan(u_exp_4)]
u_exp_6 = u_exp['x+6']
u_exp_6_cleaned = u_exp_6[~np.isnan(u_exp_6)]
u_exp_10 = u_exp['x+10']
u_exp_10_cleaned = u_exp_10[~np.isnan(u_exp_10)]

u_moose_1 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_1_sampler_0002.csv']
u_moose_4 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_4_sampler_0002.csv']
u_moose_6 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_6_sampler_0002.csv']
u_moose_10 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_10_sampler_0002.csv']

if plot_current_results:
    sim_u_moose_1 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv']
    sim_u_moose_4 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_4_sampler_0002.csv']
    sim_u_moose_6 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_6_sampler_0002.csv']
    sim_u_moose_10 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_10_sampler_0002.csv']


### Error calculations based on (ERCOFTAC - MOOSE FV K-epsilon Reference)
error_magnitude_1 = abs((1 + err_max) * (u_exp_1_cleaned - u_moose_1))
error_magnitude_4 = abs((1 + err_max) * (u_exp_4_cleaned - u_moose_4))
error_magnitude_6 = abs((1 + err_max) * (u_exp_6_cleaned - u_moose_6))
error_magnitude_10 = abs((1 + err_max) * (u_exp_10_cleaned - u_moose_10))

min_error_1 = u_exp_1 - error_magnitude_1
max_error_1 = u_exp_1 + error_magnitude_1
min_error_4 = u_exp_4 - error_magnitude_4
max_error_4 = u_exp_4 + error_magnitude_4
min_error_6 = u_exp_6 - error_magnitude_6
max_error_6 = u_exp_6 + error_magnitude_6
min_error_10 = u_exp_10 - error_magnitude_10
max_error_10 = u_exp_10 + error_magnitude_10


# x/H = 1 and x/H = 4
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
fig.suptitle("Vertical X-Velocity Profiles with Error Bars", fontsize="14")

ax1.errorbar(u_exp_1_cleaned, y_grid, xerr=error_magnitude_1, ecolor='r', capsize = 3, fmt='k.', barsabove=False)
ax1.plot(u_exp_1_cleaned, y_grid, marker='o', markersize=4, label='ERCOFTAC')
if plot_reference_results:
    ax1.plot(u_moose_1, y_grid, linestyle='--', color='g', marker='o', markersize=4, label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax1.plot(sim_u_moose_1, y_grid, linestyle='--', color='b', marker='o', markersize=4, label='MOOSE-FV K-epsilon current')
ax1.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax1.legend()
ax1.grid(True)
ax1.set_title("x/H = 1")
ax1.set_ylim([-0.1, 4.5])

ax2.errorbar(u_exp_4_cleaned, y_grid, xerr=error_magnitude_4, ecolor='r', capsize = 3, fmt='k.', barsabove=False)
ax2.plot(u_exp_4_cleaned, y_grid, marker='o', markersize=4, label='ERCOFTAC')
if plot_reference_results:
    ax2.plot(u_moose_4, y_grid, linestyle='--', color='g', marker='o', markersize=4, label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax2.plot(sim_u_moose_4, y_grid, linestyle='--', color='b', marker='o', markersize=4, label='MOOSE-FV K-epsilon current')
ax2.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax2.legend()
ax2.grid(True)
ax2.set_title("x/H = 4")
ax2.set_ylim([-0.1, 4.5])

plt.savefig('plots_u_profiles_error1_main.png', dpi=300)


# x/H = 6 and x/H = 10
fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
fig.suptitle("Vertical X-Velocity Profiles with Error Bars", fontsize="14")

ax3.errorbar(u_exp_6_cleaned, y_grid, xerr=error_magnitude_6, ecolor='r', capsize = 3, fmt='k.', barsabove=False)
ax3.plot(u_exp_6_cleaned, y_grid, marker='o', markersize=4, label='ERCOFTAC')
if plot_reference_results:
    ax3.plot(u_moose_6, y_grid, linestyle='--', color='g', marker='o', markersize=4, label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax3.plot(sim_u_moose_6, y_grid, linestyle='--', color='b', marker='o', markersize=4, label='MOOSE-FV K-epsilon current')
ax3.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax3.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax3.legend()
ax3.grid(True)
ax3.set_title("x/H = 6")
ax3.set_ylim([-0.1, 4.5])

ax4.errorbar(u_exp_10_cleaned, y_grid, xerr=error_magnitude_10, ecolor='r', capsize = 3, fmt='k.', barsabove=False)
ax4.plot(u_exp_10_cleaned, y_grid, marker='o', markersize=4, label='ERCOFTAC')
if plot_reference_results:
    ax4.plot(u_moose_10, y_grid, linestyle='--', color='g', marker='o', markersize=4, label='MOOSE-FV K-epsilon' + ref_suffix)
if plot_current_results:
    ax4.plot(sim_u_moose_10, y_grid, linestyle='--', color='b', marker='o', markersize=4, label='MOOSE-FV K-epsilon current')
ax4.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax4.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax4.legend()
ax4.grid(True)
ax4.set_title("x/H = 10")
ax4.set_ylim([-0.1, 4.5])

plt.savefig('plots_u_profiles_error2_main.png', dpi=300)
