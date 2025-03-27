#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

os.chdir(os.path.dirname(os.path.realpath(__file__)))

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12

csv_file = '../../../../../../../../validation/free_flow/isothermal/vortex_shedding/cylinder/gold/strouhal.csv'
data = pd.read_csv(csv_file)
value = data.iloc[0, 0]

# Define the acceptable range
min_range = 0.295
max_range = 0.305

fig, ax = plt.subplots()
ax.plot(1, value, 'X', color='red', markersize=10, label='Current value')
ax.errorbar(1, value, yerr=[[value - min_range], [max_range - value]], fmt='none', color='k', capsize=5, capthick=2, label='Acceptable Range')

ax.grid(True, which='major', linestyle='-', linewidth='0.5', color='grey')
ax.get_xaxis().set_visible(False)
ax.set_ylim(min_range - 0.005, max_range + 0.005)
ax.set_ylabel('Strouhal Number')

ax.legend()
plt.tight_layout()

# Show the plot
plt.savefig("strouhal.png")
