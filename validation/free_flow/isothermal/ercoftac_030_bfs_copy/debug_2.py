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

moose_cp_x = x/H
moose_cp = (pressure / cp_factor) + 0.125
ercoftac_cp = ercoftac_csv['cp']


### Interpolate MOOSE onto ERCOFTAC
moose_cp_interp = np.interp(ercoftac_csv['x/h'], moose_cp_x, moose_cp)

### Root Mean Square Error
difference_3 = moose_cp_interp - ercoftac_cp
square_3 = difference_3 ** 2
summation_3 = sum(square_3)
divide_3 = summation_3 / len(moose_cp_x)
rmse_value_cp_3 = divide_3 ** 0.5 # Calculated RMSE

# print("Pressure Coefficient RMSE: ", rmse_value_cp_3)

range_cp = max(ercoftac_cp) - min(ercoftac_cp)
normalized_rmse_cf = round((rmse_value_cp_3 / range_cp) * 100, 3)
print("Pressure Coefficient NRMSE: ", normalized_rmse_cf, " %")







error_magnitude = (abs( 1.02 * (ercoftac_cp - moose_cp_interp) ) )


# Error Ranges
min_error = ercoftac_cp - error_magnitude
max_error = ercoftac_cp + error_magnitude
# print(error_magnitude)
# print(error_bar_magnitude)

yes_no = []

for i in range(len(ercoftac_cp)):
    cp_data = moose_cp_interp[i]
    min_error_2 = min_error[i]
    max_error_2 = max_error[i]

    if (cp_data >= min_error_2) and (cp_data <= max_error_2):
        yes_no.append(0)
    else:
        yes_no.append(1)

print(yes_no)

value = 0

if sum(yes_no) == 0:
    value = 2
else:
    value = 10

print(value)









fig, (ax2) = plt.subplots(1, 1, figsize=(10, 5), constrained_layout=True)
ax2.errorbar(ercoftac_csv['x/h'], ercoftac_cp, yerr=error_magnitude, ecolor='r', capsize = 3, fmt='k--', barsabove=False)
ax2.plot(ercoftac_csv['x/h'], ercoftac_cp, 'k.', markersize=4, label='Exp')
ax2.plot(ercoftac_csv['x/h'], moose_cp_interp, linestyle='--', marker='o', markersize=4, label='interp')
ax2.set_xlabel(r'$\mathrm{x/h}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{c_p}$', fontsize=14)
ax2.set_xlim([-5,22])
ax2.legend()
ax2.grid(True)
ax2.set_title("Interp. MOOSE onto ERCOFTAC")
plt.show()