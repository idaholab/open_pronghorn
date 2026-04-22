from TestHarness.validation import ValidationCase

import os
import subprocess
import sys
import numpy as np
import pandas as pd


def csv_matches_gold(current, gold, rel_err, abs_zero):
    if list(current.columns) != list(gold.columns) or len(current) != len(gold):
        return False

    for column in current.columns:
        current_values = current[column]
        gold_values = gold[column]
        current_numeric = pd.to_numeric(current_values, errors="coerce")
        gold_numeric = pd.to_numeric(gold_values, errors="coerce")

        if current_numeric.notna().all() and gold_numeric.notna().all():
            current_array = current_numeric.to_numpy(dtype=float)
            gold_array = gold_numeric.to_numpy(dtype=float)
            current_array[np.abs(current_array) < abs_zero] = 0.0
            gold_array[np.abs(gold_array) < abs_zero] = 0.0
            if not np.allclose(
                current_array,
                gold_array,
                rtol=rel_err,
                atol=abs_zero,
                equal_nan=False,
            ):
                return False
        elif not current_values.equals(gold_values):
            return False

    return True


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

        rel_err = float(self.getParam("rel_err"))
        abs_zero = float(self.getParam("abs_zero"))
        if csv_matches_gold(
            current_centerline, gold_centerline, rel_err, abs_zero
        ) and csv_matches_gold(current_halfwidth, gold_halfwidth, rel_err, abs_zero):
            repo_root = os.path.abspath(
                os.path.join(os.path.dirname(__file__), "../../../..")
            )
            plot_script = os.path.join(
                repo_root,
                "doc/content/media/validation/jets/single-stream-jet/plot_results.py",
            )
            subprocess.run([sys.executable, plot_script], check=True)

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

        self.min_error_centerline = (
            1.0 + self.lower_bound
        ) * self.reference_sim_error_centerline
        self.max_error_centerline = (
            1.0 + self.upper_bound
        ) * self.reference_sim_error_centerline

        self.min_error_halfwidth = (
            1.0 + self.lower_bound
        ) * self.reference_sim_error_halfwidth
        self.max_error_halfwidth = (
            1.0 + self.upper_bound
        ) * self.reference_sim_error_halfwidth

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam("validation_lower_bound", "Lower bound")
        params.addRequiredParam("validation_upper_bound", "Upper bound")
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
