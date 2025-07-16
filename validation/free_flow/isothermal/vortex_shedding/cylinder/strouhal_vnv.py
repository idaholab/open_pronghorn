from TestHarness.validation import ValidationCase
import numpy as np
import pandas as pd
import os

class TestCase(ValidationCase):
    def initialize(self):

        self.gold = pd.read_csv('gold/strouhal.csv').iloc[0, 0]

        self.lower_bound = float(self.getParam('validation_lower_bound'))
        self.upper_bound = float(self.getParam('validation_upper_bound'))

        self.regression_rel_err = float(self.getParam('regression_rel_err'))


        time_series = pd.read_csv("flow_csv.csv")
        signal = time_series["lift_coeff"]
        t = time_series["time"]

        d_cylinder = 0.1
        u_bulk = 1.0

        fft_result = np.fft.fft(signal)

        frequencies = np.fft.fftfreq(len(fft_result), d=0.001)
        magnitude = np.abs(fft_result)

        peak_index = np.argmax(magnitude[:len(magnitude)//2])
        vortex_shedding_frequency = frequencies[peak_index]

        self.value = (d_cylinder * vortex_shedding_frequency / u_bulk)

        st_df = pd.DataFrame([self.value], columns=["Strouhal"])
        st_df.to_csv("strouhal.csv", index=False)

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam('validation_lower_bound', 'The lower bound for the data')
        params.addRequiredParam('validation_upper_bound', 'The upper bound for the data')
        params.addRequiredParam('regression_rel_err', 'The lower bound for the regression test')
        return params

    def testValidation(self):
        self.addScalarData("strouhal", self.value, "Strouhal number", None,
                           bounds=(self.lower_bound,self.upper_bound))

    def testVerification(self):
        self.addScalarData('strouhal_regression', self.value, 'Strouhal number', None,
                           validation=False,
                           bounds=(self.gold*(1.0-self.regression_rel_err),self.gold*(1.0+self.regression_rel_err)))
