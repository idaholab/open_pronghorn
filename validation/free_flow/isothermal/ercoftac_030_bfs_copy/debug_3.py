import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ### Load .csv files for MOOSE and ERCOFTAC, respectively
# moose_inlet_csv = pd.read_csv('bfs_input_csv_inlet_sampler_0002.csv')
# moose_outlet_csv = pd.read_csv('bfs_input_csv_outlet_sampler_0002.csv')
# ercoftac_csv = pd.read_csv('csv/cf.csv')

# # self.lower_bound = float(self.getParam('validation_lower_bound'))
# # self.upper_bound = float(self.getParam('validation_upper_bound'))

# ### Concatenate the MOOSE data at the inlet and outlet
# x = np.concatenate([moose_inlet_csv['x'], moose_outlet_csv['x']])
# mu_t = np.concatenate([moose_inlet_csv['mu_t'], moose_outlet_csv['mu_t']]) / 6.0
# distance = np.concatenate([moose_inlet_csv['distance'], moose_outlet_csv['distance']])
# vel_x = np.concatenate([moose_inlet_csv['vel_x'], moose_outlet_csv['vel_x']])

# H = 0.0127
# cf_factor = 0.5*1.18415*47.**2 # 1/2 rho U_ref^2 

# moose_cf_x = x/H
# moose_cf = (mu_t*vel_x/distance)/cf_factor
# ercoftac_cf = ercoftac_csv['cf']


# ### Interpolate MOOSE onto ERCOFTAC
# moose_cf_interp = np.interp(ercoftac_csv['x/h'], moose_cf_x, moose_cf)
# # self.moose_cf = moose_cf_interp


# error_magnitude = abs(1.02 * (ercoftac_cf - moose_cf_interp)) 

# min_error = ercoftac_cf - error_magnitude
# max_error = ercoftac_cf + error_magnitude
# # print(error_magnitude)
# # print(error_bar_magnitude)


# # # Save the output to CSV
# # df = pd.DataFrame([self.value], columns=["RMSE_Ux"])
# # df.to_csv("RMSE_Ux_395.csv", index=False)
# yes_no = []

# for i in range(len(ercoftac_cf)):
#     cf_data = moose_cf_interp[i]
#     min_error_2 = min_error[i]
#     max_error_2 = max_error[i]

#     if (cf_data >= min_error_2) and (cf_data <= max_error_2):
#         yes_no.append(0)
#     else:
#         yes_no.append(1)

# value = 0

# if sum(yes_no) == 0:
#     value = 2
# else:
#     value = 10

# print(value)





### Read .csv files and retrieve important data
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

# Extract necessary data
u_exp_1 = u_exp['x+1']
u_exp_1_cleaned = u_exp_1[~np.isnan(u_exp_1)]
u_exp_4 = u_exp['x+4']
u_exp_4_cleaned = u_exp_4[~np.isnan(u_exp_4)]
u_exp_6 = u_exp['x+6']
u_exp_6_cleaned = u_exp_6[~np.isnan(u_exp_6)]
u_exp_10 = u_exp['x+10']
u_exp_10_cleaned = u_exp_10[~np.isnan(u_exp_10)]
u_moose_1 = interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv']
u_moose_4 = interpolated_results['bfs_input_csv_vel_x_xoH_4_sampler_0002.csv']
u_moose_6 = interpolated_results['bfs_input_csv_vel_x_xoH_6_sampler_0002.csv']
u_moose_10 = interpolated_results['bfs_input_csv_vel_x_xoH_10_sampler_0002.csv']


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


yes_no_1 = [] # Initialize an empty array
yes_no_4 = []
yes_no_6 = []
yes_no_10 = []

# x/H = 1
for i in range(len(y_grid)):
    moose_vel_data_1 = u_moose_1[i]
    min_error_1_ = min_error_1[i]
    max_error_1_ = max_error_1[i]

    if (moose_vel_data_1 >= min_error_1_) and (moose_vel_data_1 <= max_error_1_):
        yes_no_1.append(0)
    else:
        yes_no_1.append(1)

# x/H = 4
for i in range(len(y_grid)):
    moose_vel_data_4 = u_moose_4[i]
    min_error_4_ = min_error_4[i]
    max_error_4_ = max_error_4[i]

    if (moose_vel_data_4 >= min_error_4_) and (moose_vel_data_4 <= max_error_4_):
        yes_no_4.append(0)
    else:
        yes_no_4.append(1)

# x/H = 6
for i in range(len(y_grid)):
    moose_vel_data_6 = u_moose_6[i]
    min_error_6_ = min_error_6[i]
    max_error_6_ = max_error_6[i]

    if (moose_vel_data_6 >= min_error_6_) and (moose_vel_data_6 <= max_error_6_):
        yes_no_6.append(0)
    else:
        yes_no_6.append(1)

# x/H = 10
for i in range(len(y_grid)):
    moose_vel_data_10 = u_moose_10[i]
    min_error_10_ = min_error_10[i]
    max_error_10_ = max_error_10[i]

    if (moose_vel_data_10 >= min_error_10_) and (moose_vel_data_10 <= max_error_10_):
        yes_no_10.append(0)
    else:
        yes_no_10.append(1)


result_array = yes_no_1 + yes_no_4 + yes_no_6 + yes_no_10

value = 0

if sum(result_array) == 0:
    value = 2
else:
    value = 10


# print(yes_no_1)
# print(yes_no_4)
# print(yes_no_6)
# print(yes_no_10)
# print(value)



### Plotting graphs with error bars

# x/H = 1 and x/H = 4
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax1.errorbar(u_exp_1_cleaned, y_grid, xerr=error_magnitude_1, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax1.plot(u_exp_1_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax1.plot(u_moose_1, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax1.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax1.legend()
ax1.grid(True)
ax1.set_title("x/H = 1")

ax2.errorbar(u_exp_4_cleaned, y_grid, xerr=error_magnitude_4, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax2.plot(u_exp_4_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax2.plot(u_moose_4, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax2.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax2.legend()
ax2.grid(True)
ax2.set_title("x/H = 4")


# x/H = 6 and x/H = 10
fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax3.errorbar(u_exp_6_cleaned, y_grid, xerr=error_magnitude_1, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax3.plot(u_exp_6_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax3.plot(u_moose_6, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax3.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax3.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax3.legend()
ax3.grid(True)
ax3.set_title("x/H = 1")

ax4.errorbar(u_exp_10_cleaned, y_grid, xerr=error_magnitude_4, ecolor='r', capsize = 3, fmt='k-', barsabove=False)
ax4.plot(u_exp_10_cleaned, y_grid, 'k.', markersize=4, label='Exp')
ax4.plot(u_moose_10, y_grid, linestyle='--', marker='o', markersize=4, label='interp')
ax4.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax4.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax4.legend()
ax4.grid(True)
ax4.set_title("x/H = 4")
plt.show()