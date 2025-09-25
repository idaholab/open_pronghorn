from TestHarness.validation import ValidationCase

import pandas as pd
import numpy as np
import os

class TestCase(ValidationCase):
    def initialize(self):

        ### Load .csv files for MOOSE and ERCOFTAC, respectively
        moose_csv = pd.read_csv('Ret590_MOOSE_LSFV_k-epsilon_test_csv_0002.csv')
        ercoftac_csv = pd.read_csv('gold/Ret590.csv')

        ### Extract necessary data from the .csv files
        # NOTE: ERCOFTAC data ordered with ascending U-mean and MOOSE is ordered with descending U-mean.
        moose_vel_x = moose_csv['vel_x'].to_numpy()
        moose_yplus = moose_csv['yplus'].iloc[-1]
        print('\n Current MOOSE Y+ value: ', moose_yplus)

        # Because the ERCOFTAC results has more data points than the MOOSE results...
        # Interpolate the fine ERCOFTAC data onto the MOOSE-measurements
        ercoftac_U_mean = np.interp(1.0 - moose_csv['y'], ercoftac_csv['y'], ercoftac_csv['Umean'])

        ### Compute friction velocity u_tau (rho = 1)
        nu = 2. / 22250
        yd = 0.08
        u_tau_moose = (moose_yplus * nu) / yd

        ### Normalize the velocities from MOOSE by dividing by friction velocity u_tau
        moose_vel_x_norm = moose_vel_x / u_tau_moose

        ### Save coordinates and result
        self.y = moose_csv['y']
        self.velocity = moose_vel_x_norm

        ### Retrieve the tolerances set in the input
        self.lower_bound = float(self.getParam('validation_lower_bound'))
        self.upper_bound = float(self.getParam('validation_upper_bound'))

        ### Get the error of the MOOSE reference results
        moose_lsfv_ref_csv = pd.read_csv('gold/Ret590_MOOSE_LSFV_k-epsilon_test_csv_0002.csv')
        moose_lsfv_ref_yp = moose_lsfv_ref_csv['yplus'].iloc[-1]
        moose_lsfv_ref_vel = moose_lsfv_ref_csv['vel_x'].to_numpy()
        print('\n Reference MOOSE Y+ value: ', moose_lsfv_ref_yp)

        moose_lsfv_vel_norm = moose_lsfv_ref_vel / (moose_lsfv_ref_yp * nu) * yd

        # Compute the errors, wrt to reference MOOSE solution and ercoftac
        self.current_sim_error = np.abs(moose_vel_x_norm - ercoftac_U_mean)
        self.reference_sim_error = np.abs(moose_lsfv_vel_norm - ercoftac_U_mean)
        self.sim_min_error = (1.0 + self.lower_bound) * self.reference_sim_error
        self.sim_max_error = (1.0 + self.upper_bound) * self.reference_sim_error
        self.yplus_rel_diff = (moose_yplus - moose_lsfv_ref_yp) / moose_lsfv_ref_yp


    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam('validation_lower_bound', 'The lower bound for the RMS on the velocity profile')
        params.addRequiredParam('validation_upper_bound', 'The upper bound for the RMS on the velocity profile')
        return params


    def testValidation(self):
        self.addVectorData('Ux_diff',
                    (self.y, 'Wall distance', '-'),
                    (self.current_sim_error, 'Error in normalized X-velocity', '-'),
                    bounds=((self.sim_min_error, self.sim_max_error)))
        self.addVectorData('Ux', (self.y, 'Wall distance', '-'), (self.velocity, 'Normalized X-velocity', '-'))
        self.addScalarData(key="y_plus", value=self.yplus_rel_diff, description="Relative difference on outlet wall adjacent cell Y plus", units='-',
                           bounds=(self.lower_bound, self.upper_bound))
