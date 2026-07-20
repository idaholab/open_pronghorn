from TestHarness.validation import ValidationCase

import numpy as np
import pandas as pd


class TestCase(ValidationCase):
    def initialize(self):
        D = 0.01

        current_centerline = (
            pd.read_csv("jet_RTL10000_out_centerline_0002.csv")
            .sort_values("z")
            .reset_index(drop=True)
        )

        current_halfwidth = (
            pd.read_csv("combined_halfwidths.csv")
            .sort_values("z")
            .reset_index(drop=True)
        )

        gold_centerline = (
            pd.read_csv("gold/jet_RTL10000_out_centerline_0002.csv")
            .sort_values("z")
            .reset_index(drop=True)
        )

        gold_halfwidth = (
            pd.read_csv("gold/combined_halfwidths.csv")
            .sort_values("z")
            .reset_index(drop=True)
        )

        exp_centerline = pd.read_csv(
            "reference/exp_centerline.csv",
            header=None,
            names=["zD", "Vin_over_Vc"],
        )

        exp_halfwidth = pd.read_csv(
            "reference/exp_halfwidth.csv",
            header=None,
            names=["zD", "rhalf_over_D"],
        )

        # Current simulation coordinates
        self.zD_centerline = current_centerline["z"].to_numpy() / D
        self.zD_halfwidth = current_halfwidth["z"].to_numpy() / D

        # Normalize centerline by first sampled point
        current_centerline_norm = (
            current_centerline["vel_z"].iloc[0] / current_centerline["vel_z"].to_numpy()
        )
        gold_centerline_norm = (
            gold_centerline["vel_z"].iloc[0] / gold_centerline["vel_z"].to_numpy()
        )

        # Normalize half-width by nozzle diameter
        current_halfwidth_norm = current_halfwidth["y_half"].to_numpy() / D
        gold_halfwidth_norm = gold_halfwidth["y_half"].to_numpy() / D

        # Interpolate experimental/reference curves onto simulation axial locations
        ref_centerline = np.interp(
            self.zD_centerline,
            exp_centerline["zD"].to_numpy(),
            exp_centerline["Vin_over_Vc"].to_numpy(),
        )

        ref_halfwidth = np.interp(
            self.zD_halfwidth,
            exp_halfwidth["zD"].to_numpy(),
            exp_halfwidth["rhalf_over_D"].to_numpy(),
        )

        self.lower_bound = float(self.getParam("validation_lower_bound"))
        self.upper_bound = float(self.getParam("validation_upper_bound"))
        validation_abs_tol = float(self.getParam("validation_abs_tol"))

        self.reference_sim_error_centerline = np.abs(
            ref_centerline - gold_centerline_norm
        )
        self.current_sim_error_centerline = np.abs(
            ref_centerline - current_centerline_norm
        )

        self.reference_sim_error_halfwidth = np.abs(ref_halfwidth - gold_halfwidth_norm)
        self.current_sim_error_halfwidth = np.abs(
            ref_halfwidth - current_halfwidth_norm
        )

        self.min_error_centerline = np.maximum(
            0.0,
            (1.0 + self.lower_bound) * self.reference_sim_error_centerline
            - validation_abs_tol,
        )
        self.max_error_centerline = (
            1.0 + self.upper_bound
        ) * self.reference_sim_error_centerline + validation_abs_tol

        self.min_error_halfwidth = np.maximum(
            0.0,
            (1.0 + self.lower_bound) * self.reference_sim_error_halfwidth
            - validation_abs_tol,
        )
        self.max_error_halfwidth = (
            1.0 + self.upper_bound
        ) * self.reference_sim_error_halfwidth + validation_abs_tol

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam("validation_lower_bound", "Lower bound")
        params.addRequiredParam("validation_upper_bound", "Upper bound")
        params.addRequiredParam(
            "validation_abs_tol",
            "Absolute tolerance for validation error bounds near zero reference error",
        )
        return params

    def testValidation(self):
        self.addVectorData(
            "centerline_diff",
            (self.zD_centerline, "z/D", "-"),
            (
                self.current_sim_error_centerline,
                "Error in normalized centerline velocity",
                "-",
            ),
            bounds=(self.min_error_centerline, self.max_error_centerline),
        )

        self.addVectorData(
            "halfwidth_diff",
            (self.zD_halfwidth, "z/D", "-"),
            (
                self.current_sim_error_halfwidth,
                "Error in normalized jet half-width",
                "-",
            ),
            bounds=(self.min_error_halfwidth, self.max_error_halfwidth),
        )
