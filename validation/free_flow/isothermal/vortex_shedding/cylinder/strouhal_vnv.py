from TestHarness.validation import ValidationCase
import numpy as np
import pandas as pd


class TestCase(ValidationCase):
    def initialize(self):

        self.gold = pd.read_csv("gold/strouhal.csv").iloc[0, 0]

        self.lower_bound = float(self.getParam("validation_lower_bound"))
        self.upper_bound = float(self.getParam("validation_upper_bound"))
        self.drag_lower_bound = float(self.getParam("drag_validation_lower_bound"))
        self.drag_upper_bound = float(self.getParam("drag_validation_upper_bound"))
        self.lift_lower_bound = float(self.getParam("lift_validation_lower_bound"))
        self.lift_upper_bound = float(self.getParam("lift_validation_upper_bound"))

        self.regression_rel_err = float(self.getParam("regression_rel_err"))

        time_series = pd.read_csv("flow_csv.csv")
        signal = time_series["lift_coeff"]
        t = time_series["time"]
        self.max_drag_coeff = time_series["drag_coeff"].max()
        self.max_lift_coeff = time_series["lift_coeff"].max()

        d_cylinder = 0.1
        u_bulk = 1.0

        fft_result = np.fft.fft(signal)

        sample_spacing = np.mean(np.diff(t))
        frequencies = np.fft.fftfreq(len(fft_result), d=sample_spacing)
        magnitude = np.abs(fft_result)

        peak_index = np.argmax(magnitude[: len(magnitude) // 2])
        vortex_shedding_frequency = frequencies[peak_index]

        self.value = d_cylinder * vortex_shedding_frequency / u_bulk

        st_df = pd.DataFrame([self.value], columns=["Strouhal"])
        st_df.to_csv("strouhal.csv", index=False)

        coefficient_df = pd.DataFrame(
            [
                {
                    "Maximum drag coefficient": self.max_drag_coeff,
                    "Maximum lift coefficient": self.max_lift_coeff,
                }
            ]
        )
        coefficient_df.to_csv("force_coefficients.csv", index=False)

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam(
            "validation_lower_bound", "The lower bound for the data"
        )
        params.addRequiredParam(
            "validation_upper_bound", "The upper bound for the data"
        )
        params.addRequiredParam(
            "drag_validation_lower_bound",
            "The literature lower bound for the maximum drag coefficient",
        )
        params.addRequiredParam(
            "drag_validation_upper_bound",
            "The literature upper bound for the maximum drag coefficient",
        )
        params.addRequiredParam(
            "lift_validation_lower_bound",
            "The literature lower bound for the maximum lift coefficient",
        )
        params.addRequiredParam(
            "lift_validation_upper_bound",
            "The literature upper bound for the maximum lift coefficient",
        )
        params.addRequiredParam(
            "regression_rel_err", "The lower bound for the regression test"
        )
        return params

    def testValidation(self):
        self.addScalarData(
            "strouhal",
            self.value,
            "Strouhal number",
            None,
            bounds=(self.lower_bound, self.upper_bound),
        )
        self.addScalarData(
            "maximum_drag_coefficient",
            self.max_drag_coeff,
            "Maximum drag coefficient",
            None,
            bounds=(self.drag_lower_bound, self.drag_upper_bound),
        )
        self.addScalarData(
            "maximum_lift_coefficient",
            self.max_lift_coeff,
            "Maximum lift coefficient",
            None,
            bounds=(self.lift_lower_bound, self.lift_upper_bound),
        )

    def testVerification(self):
        self.addScalarData(
            "strouhal_regression",
            self.value,
            "Strouhal number",
            None,
            validation=False,
            bounds=(
                self.gold * (1.0 - self.regression_rel_err),
                self.gold * (1.0 + self.regression_rel_err),
            ),
        )
