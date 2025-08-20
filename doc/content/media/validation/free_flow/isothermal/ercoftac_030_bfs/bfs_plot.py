### ERCOFTAC Case 030: BFS with Inclined Opposite Wall
### Plot Script courtesy of Dr. Mauricio Tanoret
### Last Edited by Hailey Tran-Kieu
### Date: August 12, 2025


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


### Read MOOSE .csv files
df_in = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_inlet_sampler_0002.csv')
df_out = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_outlet_sampler_0002.csv')



############################ Pressure Coefficient Analysis ############################

# Read Benchmark .csv file
cp_exp = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/csv/cp.csv')

H = 0.0127
cp_factor = 0.5*1.18415*48.18**2 # 1/2 rho U_ref^2 

# Create subplots
fig, ax1 = plt.subplots(1, 1, figsize=(4, 4))

ax1.plot(cp_exp['x/h'], cp_exp['cp'], 'k.', label='Exp')

x = np.concatenate([df_in['x'], df_out['x']])
pressure = np.concatenate([df_in['pressure'], df_out['pressure']])
ax1.plot(x/H, pressure/cp_factor + 0.125, 'g', label='MOOSE-NS-Linear')

ax1.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax1.set_xlim([-5,22])
ax1.legend()
ax1.grid(True)

# Adjust layout
plt.title('Pressure Coefficient')
plt.tight_layout()

# plt.savefig('plots_cp.png')





############################ Wall-Skin Friction Coefficient Analysis ############################

# Read Benchmark & StarCCM+ .csv files
cf_exp = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/csv/cf.csv')

H = 0.0127
cf_factor = 0.5*1.18415*47.**2

# Create subplots
fig, ax1 = plt.subplots(1, 1, figsize=(4, 4))

ax1.plot(cf_exp['x/h'], cf_exp['cf'], 'k.', label='Exp')

x = np.concatenate([df_in['x'], df_out['x']])
mu_t = np.concatenate([df_in['mu_t'], df_out['mu_t']]) / 6.0
distance = np.concatenate([df_in['distance'], df_out['distance']])
vel_x = np.concatenate([df_in['vel_x'], df_out['vel_x']])
ax1.plot(x/H, (mu_t*vel_x/distance)/cf_factor, 'g', label='MOOSE-NS-Linear')

ax1.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{c_f}$', fontsize=14)
ax1.set_xlim([-5,22])
ax1.legend()
ax1.grid(True)

# Adjust layout
plt.title('Skin Friction Coefficient')
plt.tight_layout()

# plt.savefig('plots_cf.png')





############################ Vertical X-Velocity Profiles Analysis ############################

# Read the y_grid values
y_grid = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/csv/u+1.csv')['y/h'].to_numpy()

# Function to interpolate vel_x to y_grid
def interpolate_vel_x(y_values, vel_x_values, y_grid):
    return np.interp(y_grid, y_values, vel_x_values)

# Read and process each linear CSV
linear_files = [
    '../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_1_sampler_0002.csv',
    '../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_4_sampler_0002.csv',
    '../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_6_sampler_0002.csv',
    '../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'
]

interpolated_results = {}

for file in linear_files:
    df_linear = pd.read_csv(file)
    y_linear = (df_linear['y'].to_numpy() + 0.0127) / 0.0127
    vel_x_linear = df_linear['vel_x'].to_numpy() / 48.2
    interpolated_vel_x = interpolate_vel_x(y_linear, vel_x_linear, y_grid)
    interpolated_results[file] = interpolated_vel_x

# Read u profiles
u_exp = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/csv/u_profiles_exp.csv')

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15), constrained_layout=True)
fig.suptitle("Vertical X-Velocity Profiles", fontsize="14")

ax1.plot(u_exp['x+1'], u_exp['y/h'], 'k.', label='Exp x/h:1')
ax1.plot(interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_1_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:1')
ax1.set_title('Vertical x-Velocity Profile for x/h = 1', fontsize=14)
ax1.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax1.set_ylim([-0.1, 4.5])
ax1.legend(fontsize="10")
ax1.grid(True)

ax2.plot(u_exp['x+4'], u_exp['y/h'], 'k.', label='Exp x/h:4')
ax2.plot(interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_4_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:4')
ax2.set_title('Vertical x-Velocity Profile for x/h = 4', fontsize=14)
ax2.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax2.set_ylim([-0.1, 4.5])
ax2.legend(fontsize="10")
ax2.grid(True)

ax3.plot(u_exp['x+6'], u_exp['y/h'], 'k.', label='Exp x/h:6')
ax3.plot(interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_6_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:6')
ax3.set_title('Vertical x-Velocity Profile for x/h = 6', fontsize=14)
ax3.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax3.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax3.set_ylim([-0.1, 4.5])
ax3.legend(fontsize="10")
ax3.grid(True)

ax4.plot(u_exp['x+10'], u_exp['y/h'], 'k.', label='Exp x/h:10')
ax4.plot(interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:10')
ax4.set_title('Vertical x-Velocity Profile for x/h = 10', fontsize=14)
ax4.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax4.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax4.set_ylim([-0.1, 4.5])
ax4.legend(fontsize="10")
ax4.grid(True)


# plt.savefig('plots_u_profiles.png')




############################ Generating Error Bar Plots ############################
ercoftac_cp = cp_exp['cp']
moose_cp = (pressure / cp_factor) + 0.125

ercoftac_cf = cf_exp['cf']
moose_cf = (mu_t*vel_x/distance)/cf_factor

moose_x = x/H
ercoftac_x_cp = cp_exp['x/h']
ercoftac_x_cf = cf_exp['x/h']

### Interpolate MOOSE onto ERCOFTAC
moose_cp_interp = np.interp(cp_exp['x/h'], moose_x, moose_cp)
error_magnitude_cp = (abs( 1.02 * (ercoftac_cp - moose_cp_interp) ) )

moose_cf_interp = np.interp(cf_exp['x/h'], moose_x, moose_cf)
error_magnitude_cf = (abs( 1.02 * (ercoftac_cf - moose_cf_interp) ) )

# Error Ranges
min_error_cp = ercoftac_cp - error_magnitude_cp
max_error_cp = ercoftac_cp + error_magnitude_cp

min_error_cf = ercoftac_cf - error_magnitude_cf
max_error_cf = ercoftac_cf + error_magnitude_cf


fig, (ax5, ax6) = plt.subplots(1, 2, figsize=(15, 5), constrained_layout=True)

ax5.errorbar(ercoftac_x_cp, ercoftac_cp, yerr=error_magnitude_cp, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax5.plot(ercoftac_x_cp, ercoftac_cp, 'k.', markersize=4, label='Exp')
ax5.plot(ercoftac_x_cp, moose_cp_interp, linestyle='--', marker='o', markersize=4, label='interp')
ax5.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax5.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax5.legend()
ax5.grid(True)
ax5.set_title("Pressure Coefficient Error Bar")

ax6.errorbar(ercoftac_x_cf, ercoftac_cf, yerr=error_magnitude_cf, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax6.plot(ercoftac_x_cf, ercoftac_cf, 'k.', markersize=4, label='Exp')
ax6.plot(ercoftac_x_cf, moose_cf_interp, linestyle='--', marker='o', markersize=4, label='interp')
ax6.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax6.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax6.legend()
ax6.grid(True)
ax6.set_title("Wall-Skin Friction Coefficient Error Bar")

# plt.savefig('cp_cf_error.png')




# Extract necessary data
u_exp_1 = u_exp['x+1']
u_exp_1_cleaned = u_exp_1[~np.isnan(u_exp_1)]
u_exp_4 = u_exp['x+4']
u_exp_4_cleaned = u_exp_4[~np.isnan(u_exp_4)]
u_exp_6 = u_exp['x+6']
u_exp_6_cleaned = u_exp_6[~np.isnan(u_exp_6)]
u_exp_10 = u_exp['x+10']
u_exp_10_cleaned = u_exp_10[~np.isnan(u_exp_10)]
u_moose_1 = interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_1_sampler_0002.csv']
u_moose_4 = interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_4_sampler_0002.csv']
u_moose_6 = interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_6_sampler_0002.csv']
u_moose_10 = interpolated_results['../../../../../../../validation/free_flow/isothermal/ercoftac_030_bfs_copy/bfs_input_csv_vel_x_xoH_10_sampler_0002.csv']


### Error calculations
error_magnitude_1 = abs(1.02 * (u_exp_1_cleaned - u_moose_1)) 
error_magnitude_4 = abs(1.02 * (u_exp_4_cleaned - u_moose_4)) 
error_magnitude_6 = abs(1.02 * (u_exp_6_cleaned - u_moose_6)) 
error_magnitude_10 = abs(1.02 * (u_exp_10_cleaned - u_moose_10)) 

min_error_1 = u_exp_1 - error_magnitude_1
max_error_1 = u_exp_1 + error_magnitude_1
min_error_4 = u_exp_4 - error_magnitude_4
max_error_4 = u_exp_4 + error_magnitude_4
min_error_6 = u_exp_6 - error_magnitude_6
max_error_6 = u_exp_6 + error_magnitude_6
min_error_10 = u_exp_10 - error_magnitude_10
max_error_10 = u_exp_10 + error_magnitude_10


# x/H = 1 and x/H = 4
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
fig.suptitle("Vertical X-Velocity Profiles Error Bars", fontsize="14")

ax1.errorbar(u_exp_1_cleaned, y_grid, xerr=error_magnitude_1, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax1.plot(u_exp_1_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax1.plot(u_moose_1, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax1.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax1.legend()
ax1.grid(True)
ax1.set_title("x/H = 1")

ax2.errorbar(u_exp_4_cleaned, y_grid, xerr=error_magnitude_4, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax2.plot(u_exp_4_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax2.plot(u_moose_4, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax2.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax2.legend()
ax2.grid(True)
ax2.set_title("x/H = 4")


# x/H = 6 and x/H = 10
fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
fig.suptitle("Vertical X-Velocity Profiles Error Bars", fontsize="14")

ax3.errorbar(u_exp_6_cleaned, y_grid, xerr=error_magnitude_6, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax3.plot(u_exp_6_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax3.plot(u_moose_6, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax3.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax3.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax3.legend()
ax3.grid(True)
ax3.set_title("x/H = 6")

ax4.errorbar(u_exp_10_cleaned, y_grid, xerr=error_magnitude_10, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax4.plot(u_exp_10_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax4.plot(u_moose_10, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax4.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax4.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax4.legend()
ax4.grid(True)
ax4.set_title("x/H = 10")

# plt.savefig('u_profiles_error.png')

plt.show()