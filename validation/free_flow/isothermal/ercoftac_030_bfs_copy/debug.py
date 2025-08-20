import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



############################ Pressure Coefficient RMSE ############################

### Load .csv files for MOOSE and ERCOFTAC, respectively
moose_inlet_csv = pd.read_csv('bfs_input_csv_inlet_sampler_0002.csv')
moose_outlet_csv = pd.read_csv('bfs_input_csv_outlet_sampler_0002.csv')
ercoftac_csv = pd.read_csv('csv/cp.csv')

# self.lower_bound = float(self.getParam('validation_lower_bound'))
# self.upper_bound = float(self.getParam('validation_upper_bound'))

### Concatenate the MOOSE data at the inlet and outlet
x = np.concatenate([moose_inlet_csv['x'], moose_outlet_csv['x']])
pressure = np.concatenate([moose_inlet_csv['pressure'], moose_outlet_csv['pressure']])

H = 0.0127
cp_factor = 0.5*1.18415*48.18**2 # 1/2 rho U_ref^2      

# Because the number of data points for MOOSE and ERCOFTAC do not match...
# Interpolate the fine ERCOFTAC data onto the MOOSE-measurements
ercoftac_cp_interp = np.interp(x/H, ercoftac_csv['x/h'], ercoftac_csv['cp'])

moose_cp_x = x/H
moose_cp = (pressure / cp_factor) + 0.125

### Root Mean Square Error
difference = moose_cp - ercoftac_cp_interp
square = difference ** 2
summation = sum(square)
divide = summation / len(x)
rmse_value_cp = divide ** 0.5 # Calculated RMSE

# print("Pressure Coefficient RMSE: ", rmse_value_cp)

range_cp = max(ercoftac_csv['cp']) - min(ercoftac_csv['cp'])
normalized_rmse_cp = round((rmse_value_cp / range_cp) * 100, 3)
# print("Pressure Coefficient NRMSE: ", normalized_rmse_cp, " %")



### Interpolate MOOSE onto ERCOFTAC
moose_interp_3 = np.interp(ercoftac_csv['x/h'], moose_cp_x, moose_cp)

### Root Mean Square Error
difference_3 = moose_interp_3 - ercoftac_csv['cp']
square_3 = difference_3 ** 2
summation_3 = sum(square_3)
divide_3 = summation_3 / len(moose_cp_x)
rmse_value_cp_3 = divide_3 ** 0.5 # Calculated RMSE

print("Pressure Coefficient RMSE: ", rmse_value_cp_3)

range_cp = max(ercoftac_csv['cp']) - min(ercoftac_csv['cp'])
normalized_rmse_cf = round((rmse_value_cp_3 / range_cp) * 100, 3)
print("Pressure Coefficient NRMSE: ", normalized_rmse_cf, " %")


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
fig.suptitle("Pressure Coefficient")
ax1.plot(ercoftac_csv['x/h'], ercoftac_csv['cp'], 'k.', label='Exp')
ax1.plot(x/H, ercoftac_cp_interp, 'r.', label='Interp')
ax1.plot(x/H, moose_cp, 'b.', label='hailey')
ax1.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
# ax1.set_xlim([-5,22])
ax1.legend()
ax1.grid(True)
ax1.set_title("Interp. ERCOFTAC onto MOOSE")

ax2.plot(ercoftac_csv['x/h'], ercoftac_csv['cp'], 'k.', label='Exp')
# ax2.plot(x/H, ercoftac_cf_interp, 'r.', label='Interp')
ax2.plot(ercoftac_csv['x/h'], moose_interp_3, 'b.', label='interp')
ax2.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
# ax2.set_xlim([-5,22])
ax2.legend()
ax2.grid(True)
ax2.set_title("Interp. MOOSE onto ERCOFTAC")

############################ Wall-Skin Friction Coefficient RMSE ############################

ercoftac_cf_csv = pd.read_csv('csv/cf.csv')

### Concatenate the MOOSE data at the inlet and outlet
x = np.concatenate([moose_inlet_csv['x'], moose_outlet_csv['x']])
mu_t = np.concatenate([moose_inlet_csv['mu_t'], moose_outlet_csv['mu_t']]) / 6.0
distance = np.concatenate([moose_inlet_csv['distance'], moose_outlet_csv['distance']])
vel_x = np.concatenate([moose_inlet_csv['vel_x'], moose_outlet_csv['vel_x']])

H = 0.0127
cf_factor = 0.5*1.18415*47.**2 # 1/2 rho U_ref^2 

# Because the number of data points for MOOSE and ERCOFTAC do not match...
# Interpolate the fine ERCOFTAC data onto the MOOSE-measurements
ercoftac_cf_interp = np.interp(x/H, ercoftac_cf_csv['x/h'], ercoftac_cf_csv['cf'])

moose_cf_x = x/H
moose_cf = (mu_t*vel_x/distance)/cf_factor

### Root Mean Square Error
difference = moose_cf - ercoftac_cf_interp
square = difference ** 2
summation = sum(square)
divide = summation / len(x)
rmse_value_cf = divide ** 0.5 # Calculated RMSE

# print("Skin Friction Coefficient RMSE: ", rmse_value_cf)

range_cf = max(ercoftac_cf_csv['cf']) - min(ercoftac_cf_csv['cf'])
normalized_rmse_cf = round((rmse_value_cf / range_cf) * 100, 3)
# print("Skin Friction Coefficient NRMSE: ", normalized_rmse_cf, " %")


# Plot
fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
fig.suptitle("Wall-Skin Friction Coefficient")
ax3.plot(ercoftac_cf_csv['x/h'], ercoftac_cf_csv['cf'], 'k', label='Exp')
ax3.plot(x/H, ercoftac_cf_interp, 'r.', label='Interp')
ax3.plot(x/H, moose_cf, 'b.', label='hailey')
ax3.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax3.set_ylabel(r'$\mathrm{c_f}$', fontsize=14)
# ax3.set_xlim([-5,22])
ax3.legend()
ax3.grid(True)
ax3.set_title("Interp. ERCOFTAC onto MOOSE")


### Interpolate MOOSE onto ERCOFTAC
moose_interp_2 = np.interp(ercoftac_cf_csv['x/h'], moose_cf_x, moose_cf)

### Root Mean Square Error
difference_2 = moose_interp_2 - ercoftac_cf_csv['cf']
square_2 = difference_2 ** 2
summation_2 = sum(square_2)
divide_2 = summation_2 / len(moose_interp_2)
rmse_value_cf_2 = divide_2 ** 0.5 # Calculated RMSE

print("Skin Friction Coefficient RMSE: ", rmse_value_cf_2)

range_cf = max(ercoftac_cf_csv['cf']) - min(ercoftac_cf_csv['cf'])
normalized_rmse_cf_2 = round((rmse_value_cf_2 / range_cf) * 100, 3)
print("Skin Friction Coefficient NRMSE: ", normalized_rmse_cf_2, " %")


# Plot
ax4.plot(ercoftac_cf_csv['x/h'], ercoftac_cf_csv['cf'], 'k.', label='Exp')
# ax4.plot(x/H, ercoftac_cf_interp, 'r.', label='Interp')
ax4.plot(ercoftac_cf_csv['x/h'], moose_interp_2, 'b.', label='interp')
ax4.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax4.set_ylabel(r'$\mathrm{c_f}$', fontsize=14)
# ax4.set_xlim([-5,22])
ax4.legend()
ax4.grid(True)
ax4.set_title("Interp. MOOSE onto ERCOFTAC")
# plt.show()



############################ Interpolate Vertical X-Velocity Profiles for RMSE ############################

### Read the y_grid values of the ERCOFTAC data
y_grid = pd.read_csv('csv/u+1.csv')['y/h'].to_numpy()

### Interpolate vel_x to y_grid
def interpolate_vel_x(y_values, vel_x_values, y_grid):
    return np.interp(y_grid, y_values, vel_x_values)

# Read and process each MOOSE .csv
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
u_exp_NaN = u_exp[~np.isnan(u_exp['x+1'])] # Remove all NaN values from the dataframe



############################ Vertical X-Velocity Profile at x/H = 1 RMSE ############################

# print(u_exp_NaN)

# print(interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv'])
# print(len(interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv']))
# print(u_exp_NaN['x+1'])
# print(len(u_exp_NaN['x+1']))

### Root Mean Square Error
difference = interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv'] - u_exp_NaN['x+1']
square = difference ** 2
summation = sum(square)
divide = summation / len(u_exp_NaN['x+1'])
rmse_value_xh1 = divide ** 0.5 # Calculated RMSE

print("x/H = 1 RMSE: ", rmse_value_xh1)

range_xh1 = max(u_exp_NaN['x+1']) - min(u_exp_NaN['x+1'])
normalized_rmse_xh1 = round((rmse_value_xh1 / range_xh1) * 100, 3) 
print("x/H = 1 NRMSE: ", normalized_rmse_xh1, " %")
# print(range_xh1)



############################ Vertical X-Velocity Profile at x/H = 4 RMSE ############################

### Root Mean Square Error
difference = interpolated_results['bfs_input_csv_vel_x_xoH_4_sampler_0002.csv'] - u_exp_NaN['x+4']
square = difference ** 2
summation = sum(square)
divide = summation / len(u_exp_NaN['x+1'])
rmse_value_xh4 = divide ** 0.5 # Calculated RMSE

print("x/H = 4 RMSE: ", rmse_value_xh4)
range_xh4 = max(u_exp_NaN['x+4']) - min(u_exp_NaN['x+4'])
normalized_rmse_xh4 = round((rmse_value_xh4 / range_xh4) * 100, 3) 
print("x/H = 4 NRMSE: ", normalized_rmse_xh4, " %")


############################ Vertical X-Velocity Profile at x/H = 6 RMSE ############################

### Root Mean Square Error
difference = interpolated_results['bfs_input_csv_vel_x_xoH_6_sampler_0002.csv'] - u_exp_NaN['x+6']
square = difference ** 2
summation = sum(square)
divide = summation / len(u_exp_NaN['x+1'])
rmse_value_xh6 = divide ** 0.5 # Calculated RMSE

print("x/H = 6 RMSE: ", rmse_value_xh6)
range_xh6 = max(u_exp_NaN['x+6']) - min(u_exp_NaN['x+6'])
normalized_rmse_xh6 = round((rmse_value_xh6 / range_xh6) * 100, 3) 
print("x/H = 6 NRMSE: ", normalized_rmse_xh6, " %")
# print(range_xh6)


############################ Vertical X-Velocity Profile at x/H = 10 RMSE ############################

### Root Mean Square Error
difference = interpolated_results['bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'] - u_exp_NaN['x+10']
square = difference ** 2
summation = sum(square)
divide = summation / len(u_exp_NaN['x+1'])
rmse_value_xh10 = divide ** 0.5 # Calculated RMSE

print("x/H = 10 RMSE: ", rmse_value_xh10)
range_xh10 = max(u_exp_NaN['x+10']) - min(u_exp_NaN['x+10'])
normalized_rmse_xh10 = round((rmse_value_xh10 / range_xh10) * 100, 3) 
print("x/H = 10 NRMSE: ", normalized_rmse_xh10, " %")
# print(range_xh10)

plt.show()