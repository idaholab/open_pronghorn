from TestHarness.validation import ValidationCase

import pandas as pd
import numpy as np
import os

class TestCase(ValidationCase):
    def initialize(self):

        # Authorized relative increase in the error
        err_max = 0.01

        ############################ Setting Up Reference Data ############################

        ### REFERENCE FILES: Please do not modify, unless updating reference data
        ### Load .csv files for MOOSE and ERCOFTAC, respectively
        moose_inlet_csv = pd.read_csv('reference_csv/reference_bfs_input_csv_inlet_sampler_0002.csv')
        moose_outlet_csv = pd.read_csv('reference_csv/reference_bfs_input_csv_outlet_sampler_0002.csv')

        ercoftac_csv = pd.read_csv('reference_csv/cp.csv')
        cf_exp = pd.read_csv('reference_csv/cf.csv')


        ############################ Pressure Coefficient ############################
        ### Concatenate the MOOSE data at the inlet and outlet
        x = np.concatenate([moose_inlet_csv['x'], moose_outlet_csv['x']])
        pressure = np.concatenate([moose_inlet_csv['pressure'], moose_outlet_csv['pressure']])

        H = 0.0127
        cp_factor = 0.5*1.18415*48.18**2 # 1/2 rho U_ref^2

        moose_cp_x = x/H
        self.moose_cp_x = moose_cp_x
        moose_cp = (pressure / cp_factor) + 0.125
        ercoftac_cp = ercoftac_csv['cp']
        self.ercoftac_cp_x = ercoftac_csv['x/h']

        ### Interpolate MOOSE data onto ERCOFTAC data
        moose_cp_interp = np.interp(ercoftac_csv['x/h'], moose_cp_x, moose_cp)
        self.moose_cp = moose_cp_interp

        error_magnitude = abs((1 + err_max) * (ercoftac_cp - moose_cp_interp))

        self.min_values_cp = ercoftac_cp - error_magnitude
        self.max_values_cp = ercoftac_cp + error_magnitude

        ### Setting Up Simulation Data
        ### Modify with respective .csv files from OpenPronghorn simulations

        sim_inlet_csv = pd.read_csv('bfs_input_csv_inlet_sampler_0002.csv')
        sim_outlet_csv = pd.read_csv('bfs_input_csv_outlet_sampler_0002.csv')

        sim_x = np.concatenate([sim_inlet_csv['x'], sim_outlet_csv['x']])
        sim_x_norm = sim_x/H
        sim_pressure = np.concatenate([sim_inlet_csv['pressure'], sim_outlet_csv['pressure']])
        sim_pressure_cp = (sim_pressure / cp_factor) + 0.125

        # Interpolate pressure coefficient onto the ercoftac grid
        self.sim_cp_interp = np.interp(ercoftac_csv['x/h'], sim_x_norm, sim_pressure_cp)


        ############################ Skin Friction Coefficient ############################

        H = 0.0127
        cf_factor = 0.5*1.18415*47.**2

        # Concatenate inlet and outlet data for MOOSE .csv
        x = np.concatenate([moose_inlet_csv['x'], moose_outlet_csv['x']])
        mu_t = np.concatenate([moose_inlet_csv['mu_t'], moose_outlet_csv['mu_t']]) / 6.0
        distance = np.concatenate([moose_inlet_csv['distance'], moose_outlet_csv['distance']])
        vel_x = np.concatenate([moose_inlet_csv['vel_x'], moose_outlet_csv['vel_x']])

        sim_x = np.concatenate([sim_inlet_csv['x'], sim_outlet_csv['x']])
        sim_mu_t = np.concatenate([sim_inlet_csv['mu_t'], sim_outlet_csv['mu_t']]) / 6.0
        sim_distance = np.concatenate([sim_inlet_csv['distance'], sim_outlet_csv['distance']])
        sim_vel_x = np.concatenate([sim_inlet_csv['vel_x'], sim_outlet_csv['vel_x']])

        # cf and x/H values
        ercoftac_cf = cf_exp['cf']
        moose_cf = (mu_t*vel_x/distance)/cf_factor
        sim_moose_cf = (mu_t*vel_x/distance)/cf_factor

        moose_x = x/H
        sim_moose_x = x/H
        self.ercoftac_x_cf = cf_exp['x/h']

        ### Interpolate MOOSE onto ERCOFTAC
        self.sim_moose_cf_interp = np.interp(cf_exp['x/h'], sim_moose_x, sim_moose_cf)
        moose_cf_interp = np.interp(cf_exp['x/h'], moose_x, moose_cf)
        error_magnitude_cf = (abs( (1 + err_max) * (ercoftac_cf - moose_cf_interp) ) )

        # Error Ranges
        self.min_values_cf = ercoftac_cf - error_magnitude_cf
        self.max_values_cf = ercoftac_cf + error_magnitude_cf


        ############################ Vertical x-Velocity Profiles ############################

        # Read the y_grid values
        y_grid = pd.read_csv('reference_csv/u+1.csv')['y/h'].to_numpy()
        self.y_grid = pd.read_csv('reference_csv/u+1.csv')['y/h'].to_numpy()

        # MOOSE Reference: Function to interpolate vel_x to y_grid
        def interpolate_vel_x(y_values, vel_x_values, y_grid):
            return np.interp(y_grid, y_values, vel_x_values)

        # MOOSE Reference: Read and process each linear CSV
        linear_files = [
            'reference_csv/reference_bfs_input_csv_vel_x_xoH_1_sampler_0002.csv',
            'reference_csv/reference_bfs_input_csv_vel_x_xoH_4_sampler_0002.csv',
            'reference_csv/reference_bfs_input_csv_vel_x_xoH_6_sampler_0002.csv',
            'reference_csv/reference_bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'
        ]

        interpolated_results = {}

        for file in linear_files:
            df_linear = pd.read_csv(file)
            y_linear = (df_linear['y'].to_numpy() + 0.0127) / 0.0127 # Normalize y-values to y/H
            vel_x_linear = df_linear['vel_x'].to_numpy() / 48.2 # Normalize velocities
            interpolated_vel_x = interpolate_vel_x(y_linear, vel_x_linear, y_grid) # Interpolate MOOSE data to y/H
            interpolated_results[file] = interpolated_vel_x

        # MOOSE Current: Function to interpolate vel_x to y_grid
        def sim_interpolate_vel_x(sim_y_values, sim_vel_x_values, sim_y_grid):
            return np.interp(y_grid, sim_y_values, sim_vel_x_values)

        # Read and process each linear CSV
        sim_linear_files = [
            'bfs_input_csv_vel_x_xoH_1_sampler_0002.csv',
            'bfs_input_csv_vel_x_xoH_4_sampler_0002.csv',
            'bfs_input_csv_vel_x_xoH_6_sampler_0002.csv',
            'bfs_input_csv_vel_x_xoH_10_sampler_0002.csv'
        ]

        sim_interpolated_results = {}

        for file in sim_linear_files:
            sim_df_linear = pd.read_csv(file)
            sim_y_linear = (sim_df_linear['y'].to_numpy() + 0.0127) / 0.0127 # Normalize y-values to y/H
            sim_vel_x_linear = sim_df_linear['vel_x'].to_numpy() / 48.2 # Normalize velocities
            sim_interpolated_vel_x = interpolate_vel_x(sim_y_linear, sim_vel_x_linear, y_grid) # Interpolate MOOSE data to y/H
            sim_interpolated_results[file] = sim_interpolated_vel_x


        # Read u profiles
        u_exp = pd.read_csv('reference_csv/u_profiles_exp.csv')

        # Extract necessary data

        # ERCOFTAC Benchmark Data
        u_exp_1 = u_exp['x+1']
        u_exp_1_cleaned = u_exp_1[~np.isnan(u_exp_1)]
        u_exp_4 = u_exp['x+4']
        u_exp_4_cleaned = u_exp_4[~np.isnan(u_exp_4)]
        u_exp_6 = u_exp['x+6']
        u_exp_6_cleaned = u_exp_6[~np.isnan(u_exp_6)]
        u_exp_10 = u_exp['x+10']
        u_exp_10_cleaned = u_exp_10[~np.isnan(u_exp_10)]
        # Note: the nans are removed for the last row with invalid data from the CSV reference, downloaded as is

        # MOOSE Reference Data
        u_moose_1 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_1_sampler_0002.csv']
        u_moose_4 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_4_sampler_0002.csv']
        u_moose_6 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_6_sampler_0002.csv']
        u_moose_10 = interpolated_results['reference_csv/reference_bfs_input_csv_vel_x_xoH_10_sampler_0002.csv']

        # MOOSE Current Data
        self.sim_u_moose_1 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_1_sampler_0002.csv']
        self.sim_u_moose_4 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_4_sampler_0002.csv']
        self.sim_u_moose_6 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_6_sampler_0002.csv']
        self.sim_u_moose_10 = sim_interpolated_results['bfs_input_csv_vel_x_xoH_10_sampler_0002.csv']


        ### Error calculations based on (ERCOFTAC - MOOSE: Reference)
        error_magnitude_1 = abs((1 + err_max) * (u_exp_1_cleaned - u_moose_1))
        error_magnitude_4 = abs((1 + err_max) * (u_exp_4_cleaned - u_moose_4))
        error_magnitude_6 = abs((1 + err_max) * (u_exp_6_cleaned - u_moose_6))
        error_magnitude_10 = abs((1 + err_max) * (u_exp_10_cleaned - u_moose_10))

        # Min & max error calculations based on (ERCOFTAC - MOOSE: Reference)
        sim_min_values_1 = u_exp_1 - error_magnitude_1
        self.sim_min_values_1_cleaned = sim_min_values_1[~np.isnan(sim_min_values_1)]
        sim_max_values_1 = u_exp_1 + error_magnitude_1
        self.sim_max_values_1_cleaned = sim_max_values_1[~np.isnan(sim_max_values_1)]

        sim_min_values_4 = u_exp_4 - error_magnitude_4
        self.sim_min_values_4_cleaned = sim_min_values_4[~np.isnan(sim_min_values_4)]
        sim_max_values_4 = u_exp_4 + error_magnitude_4
        self.sim_max_values_4_cleaned = sim_max_values_4[~np.isnan(sim_max_values_4)]

        sim_min_values_6 = u_exp_6 - error_magnitude_6
        self.sim_min_values_6_cleaned = sim_min_values_6[~np.isnan(sim_min_values_6)]
        sim_max_values_6 = u_exp_6 + error_magnitude_6
        self.sim_max_values_6_cleaned = sim_max_values_6[~np.isnan(sim_max_values_6)]

        sim_min_values_10 = u_exp_10 - error_magnitude_10
        self.sim_min_values_10_cleaned = sim_min_values_10[~np.isnan(sim_min_values_10)]
        sim_max_values_10 = u_exp_10 + error_magnitude_10
        self.sim_max_values_10_cleaned = sim_max_values_10[~np.isnan(sim_max_values_10)]

    def testValidation(self):
        self.addVectorData('xh_cp',
                    (self.ercoftac_cp_x, 'Normalized distance', '-'),
                    (self.sim_cp_interp, 'Pressure coefficient', '-'),
                    bounds=((self.min_values_cp, self.max_values_cp)))

        self.addVectorData('xh_cf',
            (self.ercoftac_x_cf, 'Normalized distance', '-'),
            (self.sim_moose_cf_interp, 'Skin friction coefficient', '-'),
            bounds=((self.min_values_cf, self.max_values_cf)))

        self.addVectorData('yh_1',
            (self.y_grid, 'Normalized distance', '-'),
            (self.sim_u_moose_1, 'Vertical x-Velocity at x/H=1', '-'),
            bounds=((self.sim_min_values_1_cleaned, self.sim_max_values_1_cleaned)))

        self.addVectorData('yh_4',
            (self.y_grid, 'Normalized distance', '-'),
            (self.sim_u_moose_4, 'Vertical x-Velocity at x/H=4', '-'),
            bounds=((self.sim_min_values_4_cleaned, self.sim_max_values_4_cleaned)))

        self.addVectorData('yh_6',
            (self.y_grid, 'Normalized distance', '-'),
            (self.sim_u_moose_6, 'Vertical x-Velocity at x/H=6', '-'),
            bounds=((self.sim_min_values_6_cleaned, self.sim_max_values_6_cleaned)))

        self.addVectorData('yh_10',
            (self.y_grid, 'Normalized distance', '-'),
            (self.sim_u_moose_10, 'Vertical x-Velocity at x/H=10', '-'),
            bounds=((self.sim_min_values_10_cleaned, self.sim_max_values_10_cleaned)))
