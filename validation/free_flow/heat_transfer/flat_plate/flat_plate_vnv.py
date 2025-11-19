from TestHarness.validation import ValidationCase

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.special import erfcx


class TestCase(ValidationCase):
    def initialize(self):

        ### Load .csv files for MOOSE and ERCOFTAC, respectively
        current_interface_temp_df = pd.read_csv(
            "cht_rob-rob_out_interface_temp_0001.csv"
        ).sort_values("x")
        current_vertical_ts_df = pd.read_csv(
            "cht_rob-rob_out_y_vs_ts_0001.csv"
        ).sort_values("y")
        current_vertical_tf_df = pd.read_csv(
            "cht_rob-rob_out_y_vs_tf_0001.csv"
        ).sort_values("y")

        gold_interface_temp_df = pd.read_csv(
            "gold/cht_rob-rob_out_interface_temp_0001.csv"
        ).sort_values("x")
        gold_vertical_ts_df = pd.read_csv(
            "gold/cht_rob-rob_out_y_vs_ts_0001.csv"
        ).sort_values("y")
        gold_vertical_tf_df = pd.read_csv(
            "gold/cht_rob-rob_out_y_vs_tf_0001.csv"
        ).sort_values("y")

        analytic_reference_interface = [
            self.analytic_approximate_solution_dht_erfcx(x, 1e-9)
            for x in gold_interface_temp_df["x"]
        ]
        analytic_reference_ts = [
            self.analytic_approximate_solution_dht_erfcx(0.1, y)
            for y in gold_vertical_ts_df["y"]
        ]
        analytic_reference_tf = [
            self.analytic_approximate_solution_dht_erfcx(0.1, y)
            for y in gold_vertical_tf_df["y"]
        ]

        self.x = gold_interface_temp_df["x"]
        self.y_s = current_vertical_ts_df["y"]
        self.y_f = current_vertical_tf_df["y"]

        ### Retrieve the tolerances set in the input
        self.lower_bound = float(self.getParam("validation_lower_bound"))
        self.upper_bound = float(self.getParam("validation_upper_bound"))

        # Compute the errors, wrt to reference MOOSE solution and analytic DHT solution
        self.reference_sim_error_interface = np.abs(
            analytic_reference_interface - gold_interface_temp_df["T_fluid"]
        )
        self.reference_sim_error_vertical_ts = np.abs(
            analytic_reference_ts - gold_vertical_ts_df["T_solid"]
        )
        self.reference_sim_error_vertical_tf = np.abs(
            analytic_reference_tf - gold_vertical_tf_df["T_fluid"]
        )

        self.current_sim_error_interface = np.abs(
            analytic_reference_interface - current_interface_temp_df["T_fluid"]
        )
        self.current_sim_error_vertical_ts = np.abs(
            analytic_reference_ts - current_vertical_ts_df["T_solid"]
        )
        self.current_sim_error_vertical_tf = np.abs(
            analytic_reference_tf - current_vertical_tf_df["T_fluid"]
        )

        self.min_error_interface = (
            1.0 + self.lower_bound
        ) * self.reference_sim_error_interface
        self.max_error_interface = (
            1.0 + self.upper_bound
        ) * self.reference_sim_error_interface

        self.min_error_vertical_ts = (
            1.0 + self.lower_bound
        ) * self.reference_sim_error_vertical_ts
        self.max_error_vertical_ts = (
            1.0 + self.upper_bound
        ) * self.reference_sim_error_vertical_ts

        self.min_error_vertical_tf = (
            1.0 + self.lower_bound
        ) * self.reference_sim_error_vertical_tf
        self.max_error_vertical_tf = (
            1.0 + self.upper_bound
        ) * self.reference_sim_error_vertical_tf

    def analytic_approximate_solution_dht_erfcx(self, x, y):
        # --- given properties/inputs ---
        T_b = 600.0
        T_inf = 1000.0
        lambda_s = 0.2876
        lambda_f = 0.06808
        b = 0.01
        rho_f = 0.3525
        u_inf = 12.0
        cp = 1142.6
        mu = 3.95e-5

        # --- auxiliary nondimensional groups ---
        Re_x = (u_inf * x * rho_f) / mu
        Pr = (cp * mu) / lambda_f
        u_x_bar = 0.66 * u_inf
        Re_x_bar = (x * u_x_bar * rho_f) / mu
        u_y_bar = 0.43 * u_inf * Re_x ** (-0.5)

        lambda_star = lambda_s / lambda_f
        K = lambda_star * x / b * (Pr * Re_x_bar) ** (-0.5)
        B = (u_y_bar / u_x_bar) * (Pr * Re_x_bar) ** 0.5
        alpha = lambda_star * y / b  # dimensionless wall-normal position into the fluid

        # --- scaled erfc arguments ---
        z1 = K - 0.5 * B + 0.5 * alpha / K
        z2 = 0.5 * B - 0.5 * alpha / K
        z3 = 0.5 * B + 0.5 * alpha / K

        # common exponential factor after erfcx rewrite (see derivation below)
        E = np.exp(-0.25 * (B - alpha / K) ** 2)

        # guard the K ≈ B case with the closed-form limit
        rel = abs(K - B) / max(1.0, abs(K), abs(B))
        if rel < 1e-8:
            # limit as K -> B
            z2_0 = 0.5 * B - 0.5 * alpha / B
            z3_0 = 0.5 * B + 0.5 * alpha / B
            E0 = np.exp(-0.25 * (B - alpha / B) ** 2)
            F = E0 * (
                0.5 * erfcx(z2_0)
                + 2.0 * erfcx(z3_0)
                + B * z3_0 * erfcx(z3_0)
                - B / np.sqrt(np.pi)
            )
        else:
            # stable general form
            A = (K - 0.5 * B) / (K - B)
            C = 0.5 / (1.0 - B / K)  # = 0.5*K/(K - B)
            F = E * (A * erfcx(z1) + 0.5 * erfcx(z2) - C * erfcx(z3))

        if y <= 0:
            T_w = self.analytic_approximate_solution_dht_erfcx(x, 1e-12)
            return T_w + (T_w - T_b) * y / b

        return T_b + (T_inf - T_b) * F

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam(
            "validation_lower_bound", "The lower bound for the simulation results"
        )
        params.addRequiredParam(
            "validation_upper_bound", "The upper bound for the simulation results"
        )
        return params

    def testValidation(self):
        self.addVectorData(
            "T_interface_diff",
            (self.x, "X", "m"),
            (self.current_sim_error_interface, "Error in interface temperature", "K"),
            bounds=((self.min_error_interface, self.max_error_interface)),
        )
        self.addVectorData(
            "Ts_diff",
            (self.y_s, "Y", "m"),
            (
                self.current_sim_error_vertical_ts,
                "Error in solid vertical temperature profile",
                "K",
            ),
            bounds=((self.min_error_vertical_ts, self.max_error_vertical_ts)),
        )
        self.addVectorData(
            "Tf_diff",
            (self.y_f, "Y", "m"),
            (
                self.current_sim_error_vertical_tf,
                "Error in fluid vertical temperature profile",
                "K",
            ),
            bounds=((self.min_error_vertical_tf, self.max_error_vertical_tf)),
        )
