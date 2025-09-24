import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir(os.path.dirname(os.path.realpath(__file__)))



####################################################################################################################################
### Ret 395: Comparison of axial velocity radial profiles at the outlet of the channel against DNS ###

# Load .csv files
df_395_MOOSE_linear = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_032_channel_flow/gold/Ret395_MOOSE_LSFV_k-epsilon_test_csv_0002.csv')
df_395_ercoftac = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_032_channel_flow/gold/Ret395.csv')

# Problem Parameters
nu = 2./14000. # Bulk Reynolds number
yd = 0.125 # Centroid distance from wall
rho = 1 # Density

# Compute u_tau for MOOSE Linear SIMPLE FV with k-epsilon
yplus_MOOSE_linear = df_395_MOOSE_linear['yplus'].iloc[-1] # y+ value at the outlet (obtained from the MOOSE simulation results)
u_tau_MOOSE_linear = yplus_MOOSE_linear * nu / yd

# Compute wall shear stress
tau_w_MOOSE_linear = u_tau_MOOSE_linear * u_tau_MOOSE_linear * rho
u_tau_dns = 392.5 * nu
tau_w_dns = u_tau_dns * u_tau_dns * rho

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5)) # 1, 1, figsize=(4, 4)

# Plotting
# The velocities are normalized by u_tau
ax1.plot(df_395_MOOSE_linear['y'][::2],
         df_395_MOOSE_linear['vel_x'][::2] / u_tau_MOOSE_linear, 'rx', markersize = 6, label = 'MOOSE: FV k-eps')
ax1.plot(1.0 - df_395_ercoftac['y'], df_395_ercoftac['Umean'], 'k-', label = 'DNS')

# Plot Settings
ax1.set_xlabel(r'$\mathrm{y^*}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{U^*}$', fontsize=14)
ax1.legend()
ax1.grid(True)
ax1.set_title(r'Axial Velocity Profiles for $\mathrm{Re_{\tau}~395}$')


####################################################################################################################################
### Ret 395: Comparison of the law of the wall against DNS ###

# Plotting
ax2.semilogx((1.0 - df_395_MOOSE_linear['y'][::2][10:]) * u_tau_MOOSE_linear / nu,
             df_395_MOOSE_linear['vel_x'][::2][10:] / u_tau_MOOSE_linear, 'rx', markersize = 6, label = 'MOOSE: FV k-eps')
ax2.semilogx(df_395_ercoftac['y+'], df_395_ercoftac['Umean'], 'k-', label = 'DNS')

# Plot Settings
ax2.set_xlabel(r'$\mathrm{y^+}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{U^*}$', fontsize=14)
ax2.legend()
ax2.grid(True)
ax2.set_title(r'Law of the Wall for $\mathrm{Re_{\tau}~395}$')

plt.tight_layout()
plt.savefig('7_Ret395_dual_plots')


####################################################################################################################################
### Ret 395: Comparing the wall shear stress against DNS ###

print('tau_w DNS: ',tau_w_dns )
print('tau_w MOOSE FV k-eps: ',tau_w_MOOSE_linear )
print('Error % MOOSE FV k-eps: ', 100*(tau_w_MOOSE_linear-tau_w_dns)/tau_w_dns)


####################################################################################################################################
### Ret 590: Comparison of axial velocity radial profiles at the outlet of the channel against DNS ###

# Load .csv files
df_590_ercoftac = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_032_channel_flow/gold/Ret590.csv')
df_590_MOOSE_linear = pd.read_csv('../../../../../../../validation/free_flow/isothermal/ercoftac_032_channel_flow/gold/Ret590_MOOSE_LSFV_k-epsilon_test_csv_0002.csv')

# Problem Parameters
nu = 2./22250. # Bulk Reynolds number
yd = 0.08 # Distance to wall-centroid
rho = 1 # Density

# Compute u_tau for MOOSE Linear SIMPLE FV with k-epsilon
yplus_MOOSE_linear = df_590_MOOSE_linear['yplus'].iloc[-1] # y+ obtained from the MOOSE simulation results
u_tau_MOOSE_linear = yplus_MOOSE_linear * nu / yd

# Compute wall shear stress
tau_w_MOOSE_linear = u_tau_MOOSE_linear * u_tau_MOOSE_linear
u_tau_dns = 587.2 * nu
tau_w_dns = u_tau_dns * u_tau_dns

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# Plotting
ax1.plot(df_590_MOOSE_linear['y'][::2],
         df_590_MOOSE_linear['vel_x'][::2] / u_tau_MOOSE_linear, 'rx', markersize = 6, label = 'MOOSE: FV k-eps')
ax1.plot(1.0 - df_590_ercoftac['y'], df_590_ercoftac['Umean'], 'k-', label = 'DNS')

# Plot Settings
ax1.set_xlabel(r'$\mathrm{y^*}$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{U^*}$', fontsize=14)
ax1.legend()
ax1.grid(True)
ax1.set_title(r'Axial Velocity Profile for $\mathrm{Re_{\tau}~590}$')


####################################################################################################################################
### Ret 590: Comparison of the law of the wall against DNS ###

# Plotting
ax2.semilogx((1.0 - df_590_MOOSE_linear['y'][1::2][10:]) * u_tau_MOOSE_linear / nu,
             df_590_MOOSE_linear['vel_x'][1::2][10:]/u_tau_MOOSE_linear, 'rx', markersize = 6,  label = 'MOOSE: FV k-eps')
ax2.semilogx(df_590_ercoftac['y+'], df_590_ercoftac['Umean'], 'k-', label = 'DNS')

# Plot Settings
ax2.set_xlabel(r'$\mathrm{y^+}$', fontsize=14)
ax2.set_ylabel(r'$\mathrm{U^*}$', fontsize=14)
ax2.legend()
ax2.grid(True)
ax2.set_title(r'$\mathrm{Re_{\tau}~590}$')

plt.tight_layout()
plt.savefig('8_Ret590_dual_plots')


####################################################################################################################################
### Ret 590: Comparing the wall shear stress against DNS ###

print('tau_w DNS: ',tau_w_dns )
print('tau_w MOOSE FV k-eps: ',tau_w_MOOSE_linear )
print('Error % MOOSE FV k-eps: ', 100 * (tau_w_MOOSE_linear - tau_w_dns) / tau_w_dns)
