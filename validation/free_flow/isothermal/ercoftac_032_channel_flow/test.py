import pandas as pd
import numpy as np

moose_csv = pd.read_csv('/Users/tranh/projects/open_pronghorn_other_contribs/validation/free_flow/isothermal/ercoftac_032_channel_flow/Ret395_linear_SIMPLE_k-epsilon_csv_test_csv_0002.csv')
ercoftac_csv = pd.read_csv('/Users/tranh/projects/open_pronghorn_other_contribs/validation/free_flow/isothermal/ercoftac_032_channel_flow/2_plots/Ret395.csv')

### Slicing the data
# ERCOFTAC data recorded least to greatest and MOOSE data was recorded greatest to least.
moose_vel_x = moose_csv['vel_x'][:70:1]
moose_vel_x_reversed = moose_vel_x[::-1]
ercoftac_U_mean = ercoftac_csv['Umean'][59:129:1]


# Computing friction velocity u_tau
nu = 2 / 14000
yd = 0.125
yplus_moose = 49.38567599
u_tau_moose = (yplus_moose * nu) / yd


# Normalize velocities
moose_vel_x_norm = moose_vel_x_reversed / u_tau_moose


### Renumbering the indexes so the dataframes align
df_moose_vel_x_norm_reversed = pd.DataFrame(moose_vel_x_norm, index = None)
df_moose_vel_x_norm_reversed_renumbered = df_moose_vel_x_norm_reversed.reset_index(drop = True)
df_ercoftac_U_mean = pd.DataFrame(ercoftac_U_mean, index = None)
df_ercoftac_U_mean_renumbered = df_ercoftac_U_mean.reset_index(drop = True)

difference = df_moose_vel_x_norm_reversed_renumbered['vel_x'] - df_ercoftac_U_mean_renumbered['Umean']


print(df_moose_vel_x_norm_reversed_renumbered)
print(df_ercoftac_U_mean_renumbered)
print(difference)