### ERCOFTAC Case 030: BFS with Inclined Opposite Wall
### Plot Script courtesy of Dr. Mauricio Tanoret
### Last Edited by Hailey Tran-Kieu
### Date: August 12, 2025


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


### Read MOOSE .csv files
df_in = pd.read_csv('bfs_input_csv_inlet_sampler_0002.csv')
df_out = pd.read_csv('bfs_input_csv_outlet_sampler_0002.csv')



############################ Pressure Coefficient Analysis ############################

# Read Benchmark & StarCCM+ .csv files
cp_exp = pd.read_csv('csv/cp.csv')
# cp_ph = pd.read_csv('csv/lower_wall_ph.csv')
# cp_star = pd.read_csv('csv/lower_wall_star_p.csv')
# cp_star = cp_star.sort_values(by='Region: Body 1.wall: Centroid[X] (m)')

H = 0.0127
cp_factor = 0.5*1.18415*48.18**2 # 1/2 rho U_ref^2 

# Create subplots
fig, ax1 = plt.subplots(1, 1, figsize=(4, 4))

ax1.plot(cp_exp['x/h'], cp_exp['cp'], 'k.', label='Exp')
# ax1.plot(cp_ph['Points_0']/H, cp_ph['pressure']/cp_factor +0.125, 'b', label='MOOSE-NS')
# ax1.plot(cp_star['Region: Body 1.wall: Centroid[X] (m)']/H, cp_star['Region: Body 1.wall: Pressure (Pa)']/cp_factor +0.125, 'r', label='Star-CCM')

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

plt.savefig('plots_cp.png')





############################ Wall-Skin Friction Coefficient Analysis ############################

# Read Benchmark & StarCCM+ .csv files
cp_exp = pd.read_csv('csv/cf.csv')
# cp_ph = pd.read_csv('csv/lower_wall_ph.csv')
# cp_star = pd.read_csv('csv/lower_wall_star.csv')
# cp_star = cp_star.sort_values(by='Region: Body 1.wall: Centroid[X] (m)')

H = 0.0127
cp_factor = 0.5*1.18415*47.**2
# Create subplots

fig, ax1 = plt.subplots(1, 1, figsize=(4, 4))

ax1.plot(cp_exp['x/h'], cp_exp['cf'], 'k.', label='Exp')
# ax1.plot(cp_ph['Points_0']/H, cp_ph['tau_w']/cp_factor, 'b', label='MOOSE-NS')
# ax1.plot(cp_star['Region: Body 1.wall: Centroid[X] (m)']/H, cp_star['Region: Body 1.wall: Wall Shear Stress[i] (Pa)']/cp_factor, 'r', label='Star-CCM')

x = np.concatenate([df_in['x'], df_out['x']])
mu_t = np.concatenate([df_in['mu_t'], df_out['mu_t']]) / 6.0
distance = np.concatenate([df_in['distance'], df_out['distance']])
vel_x = np.concatenate([df_in['vel_x'], df_out['vel_x']])
ax1.plot(x/H, (mu_t*vel_x/distance)/cp_factor, 'g', label='MOOSE-NS-Linear')

ax1.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{c_f}$', fontsize=14)
ax1.set_xlim([-5,22])
ax1.legend()
ax1.grid(True)

# Adjust layout
plt.title('Skin Friction Coefficient')
plt.tight_layout()

plt.savefig('plots_cf.png')





############################ Vertical X-Velocity Profiles Analysis ############################

# Read the y_grid values
y_grid = pd.read_csv('csv/u+1.csv')['y/h'].to_numpy()

# Function to interpolate vel_x to y_grid
def interpolate_vel_x(y_values, vel_x_values, y_grid):
    return np.interp(y_grid, y_values, vel_x_values)

# Read and process each linear CSV
linear_files = [
    'bfs_input_csv_vel_x_xoH_1_sampler_0002.csv',
    'bfs_input_csv_vel_x_xoH_4_sampler_0002.csv',
    'bfs_input_csv_vel_x_xoH_6_sampler_0002.csv',
    'bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'
]

interpolated_results = {}

for file in linear_files:
    df_linear = pd.read_csv(file)
    y_linear = (df_linear['y'].to_numpy() + 0.0127) / 0.0127
    vel_x_linear = df_linear['vel_x'].to_numpy() / 48.2
    interpolated_vel_x = interpolate_vel_x(y_linear, vel_x_linear, y_grid)
    interpolated_results[file] = interpolated_vel_x

# Read u profiles
u_exp = pd.read_csv('csv/u_profiles_exp.csv')
# u_ph_1 = pd.read_csv('csv/u+1_ph.csv')
# u_ph_4 = pd.read_csv('csv/u+4_ph.csv')
# u_ph_6 = pd.read_csv('csv/u+6_ph.csv')
# u_ph_10 = pd.read_csv('csv/u+10_ph.csv')
# u_star_1 = pd.read_csv('csv/u+1_star.csv').sort_values(by='Line Probe: Centroid[Y] (m)')
# u_star_4 = pd.read_csv('csv/u+4_star.csv').sort_values(by='Line Probe: Centroid[Y] (m)')
# u_star_6 = pd.read_csv('csv/u+6_star.csv').sort_values(by='Line Probe: Centroid[Y] (m)')
# u_star_10 = pd.read_csv('csv/u+10_star.csv').sort_values(by='Line Probe: Centroid[Y] (m)')

# Create subplots
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 15), constrained_layout=True)
fig.suptitle("Vertical X-Velocity Profiles", fontsize="14")

ax1.plot(u_exp['x+1'], u_exp['y/h'], 'k.', label='Exp x/h:1')
# ax1.plot(u_ph_1['vel_0']/48.2, 1 + u_ph_1['Points_1']/0.0127, 'b', label='MOOSE-NS x/h:1')
# ax1.plot(u_star_1['Line Probe: Velocity[i] (m/s)']/48.2, 1 + u_star_1['Line Probe: Centroid[Y] (m)']/0.0127, 'r', label='Star-CCM x/h:1')
ax1.plot(interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:1')
ax1.set_title('Vertical x-Velocity Profile for x/h = 1', fontsize=14)
ax1.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax1.set_ylim([-0.1, 4.5])
ax1.legend(fontsize="10")
ax1.grid(True)

ax2.plot(u_exp['x+4'], u_exp['y/h'], 'k.', label='Exp x/h:4')
# ax2.plot(u_ph_4['vel_0']/48.2, 1 + u_ph_4['Points_1']/0.0127, 'b', label='MOOSE-NS x/h:4')
# ax2.plot(u_star_4['Line Probe: Velocity[i] (m/s)']/48.2, 1 + u_star_4['Line Probe: Centroid[Y] (m)']/0.0127, 'r', label='Star-CCM x/h:4')
ax2.plot(interpolated_results['bfs_input_csv_vel_x_xoH_4_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:4')
ax2.set_title('Vertical x-Velocity Profile for x/h = 4', fontsize=14)
ax2.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax2.set_ylim([-0.1, 4.5])
ax2.legend(fontsize="10")
ax2.grid(True)

ax3.plot(u_exp['x+6'], u_exp['y/h'], 'k.', label='Exp x/h:6')
# ax3.plot(u_ph_6['vel_0']/48.2, 1 + u_ph_6['Points_1']/0.0127, 'b', label='MOOSE-NS x/h:6')
# ax3.plot(u_star_6['Line Probe: Velocity[i] (m/s)']/48.2, 1 + u_star_6['Line Probe: Centroid[Y] (m)']/0.0127, 'r', markersize=2, label='Star-CCM x/h:6')
ax3.plot(interpolated_results['bfs_input_csv_vel_x_xoH_6_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:6')
ax3.set_title('Vertical x-Velocity Profile for x/h = 6', fontsize=14)
ax3.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax3.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax3.set_ylim([-0.1, 4.5])
ax3.legend(fontsize="10")
ax3.grid(True)

ax4.plot(u_exp['x+10'], u_exp['y/h'], 'k.', label='Exp x/h:10')
# ax4.plot(u_ph_10['vel_0']/48.2, 1 + u_ph_10['Points_1']/0.0127, 'b', label='MOOSE-NS x/h:10')
# ax4.plot(u_star_10['Line Probe: Velocity[i] (m/s)']/48.2, 1 + u_star_10['Line Probe: Centroid[Y] (m)']/0.0127, 'r', markersize=2, label='Star-CCM x/h:10')
ax4.plot(interpolated_results['bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'], y_grid, 'g', label='Linear x/h:10')
ax4.set_title('Vertical x-Velocity Profile for x/h = 10', fontsize=14)
ax4.set_xlabel(r'$\mathrm{U/U_{ref}}$', fontsize=14)
ax4.set_ylabel(r'$\mathrm{y/h}$', fontsize=14)
ax4.set_ylim([-0.1, 4.5])
ax4.legend(fontsize="10")
ax4.grid(True)


# Adjust layout
plt.tight_layout()

plt.savefig('plots_u_profiles.png')