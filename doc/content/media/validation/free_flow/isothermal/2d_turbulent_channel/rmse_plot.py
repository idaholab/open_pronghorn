import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

os.chdir(os.path.dirname(os.path.realpath(__file__)))

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12

# /Users/tranh/projects/open_pronghorn_other_contribs/validation/free_flow/isothermal/vortex_shedding/cylinder/gold/strouhal.csv
csv_file_395 = '../../../../../../../validation/free_flow/isothermal/ercoftac_032_channel_flow/RMSE_Ux_395.csv'
data_395 = pd.read_csv(csv_file_395)
value_395 = data_395.iloc[0, 0]

csv_file_590 = '../../../../../../../validation/free_flow/isothermal/ercoftac_032_channel_flow/RMSE_Ux_590.csv'
data_590 = pd.read_csv(csv_file_590)
value_590 = data_590.iloc[0, 0]

# Define the acceptable range
min_range = 0.00
max_range = 0.25

fig, ax = plt.subplots()
ax.plot(1, value_395, 'X', color='red', markersize=10, label='$Re_t$ 395 Current value')
ax.plot(1, value_590, 'X', color='blue', markersize=10, label='$Re_t$ 590 Current value')
ax.errorbar(1, value_395, yerr=[[value_395 - min_range], [max_range - value_395]], fmt='none', color='k', capsize=5, capthick=2, label='Acceptable Range')

ax.grid(True, which='major', linestyle='-', linewidth='0.5', color='grey')
ax.get_xaxis().set_visible(False)
ax.set_ylim(min_range - 0.1, max_range + 0.1)
ax.set_ylabel('RMSE')

ax.legend()
plt.tight_layout()

# Show the plot
plt.savefig("rmse.png")