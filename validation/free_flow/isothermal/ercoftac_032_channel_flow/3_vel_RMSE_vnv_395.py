from TestHarness.validation import ValidationCase

import pandas as pd
import numpy as np
import os

class TestCase(ValidationCase):
    def initialize(self):

        ### Load .csv files for MOOSE and ERCOFTAC, respectively
        moose_csv = pd.read_csv('1_input/Ret395_linear_SIMPLE_k-epsilon_csv_test_csv_0002.csv')
        ercoftac_csv = pd.read_csv('2_plots/Ret395.csv')

        self.lower_bound = float(self.getParam('validation_lower_bound'))
        self.upper_bound = float(self.getParam('validation_upper_bound'))


        ### Extract necessary data from the .csv files
        # NOTE: ERCOFTAC data ordered with ascending U-mean and MOOSE is ordered with descending U-mean.
        moose_vel_x = moose_csv['vel_x'].to_numpy()
        moose_yplus = moose_csv['yplus'].iloc[-1]
        print('\n MOOSE Y+ value: ', moose_yplus)

        # Because the ERCOFTAC results has more data points than the MOOSE results...
        # Interpolate the fine ERCOFTAC data onto the MOOSE-measurements
        ercoftac_U_mean = np.interp(1.0 - moose_csv['y'], ercoftac_csv['y'], ercoftac_csv['Umean'])


        ### Compute friction velocity u_tau (rho = 1)
        nu = 2. / 14000 # Bulk Reynolds number
        yd = 0.125 # Centroid distance to wall
        u_tau_moose = (moose_yplus * nu) / yd


        ### Normalize the velocities from MOOSE by dividing by friction velocity u_tau
        moose_vel_x_norm = moose_vel_x / u_tau_moose


        ### Root Mean Square Error
        difference = moose_vel_x_norm - ercoftac_U_mean
        square = difference ** 2
        summation = sum(square)
        divide = summation / len(moose_vel_x)
        self.value = divide ** 0.5 # Calculated RMSE


        # Save the output to CSV
        df = pd.DataFrame([self.value], columns=["RMSE_Ux"])
        df.to_csv("RMSE_Ux_395.csv", index=False)


    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam('validation_lower_bound', 'The lower bound for the RMS on the velocity profile')
        params.addRequiredParam('validation_upper_bound', 'The upper bound for the RMS on the velocity profile')
        return params


    def testValidation(self):
        self.addScalarData(key="number", value=self.value, description="Root Mean Square Error", units='-',
                           bounds=(self.lower_bound,self.upper_bound))

# if __name__ == "__main__":

#     print("here")

#     a = TestCase()
#     a.initialize()
