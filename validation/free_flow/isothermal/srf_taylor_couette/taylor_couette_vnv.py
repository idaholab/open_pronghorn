from pathlib import Path

from TestHarness.validation import ValidationCase

import numpy as np
import pandas as pd


class TestCase(ValidationCase):
    """
    Validate the reconstructed absolute-frame Taylor-Couette tangential
    velocity against the analytical solution at the same FV cell centroids.

    The radial centroid row closest to y = 0 on the +x side of the annulus is
    used. Along this row,

        u_theta = vel_abs_y
        u_r     = vel_abs_x

    Pressure is also compared with the analytical radial pressure distribution
    after removing the arbitrary additive pressure constant. The pressure and
    radial-velocity comparisons are written to a CSV for inspection, while the
    automated validation criterion is applied to u_theta.
    """

    # Benchmark parameters
    RI = 0.35
    RO = 1.0
    OMEGA_I = -0.01
    OMEGA_O = 0.0
    RHO = 1.0

    # Same location used to pin pressure in the input file
    PRESSURE_REFERENCE_RADIUS = 0.75

    def initialize(self):
        current_file = self._find_centroid_csv()
        current = pd.read_csv(current_file)

        required_columns = {
            "x",
            "y",
            "vel_abs_x",
            "vel_abs_y",
            "pressure",
        }

        missing = required_columns.difference(current.columns)
        if missing:
            raise RuntimeError(
                f"{current_file} is missing required columns: {sorted(missing)}"
            )

        # ------------------------------------------------------------------
        # Select one complete radial row of FV centroids. Use the row whose
        # centroid angle is closest to the +x axis.
        # ------------------------------------------------------------------
        x_all = current["x"].to_numpy(dtype=float)
        y_all = current["y"].to_numpy(dtype=float)
        theta_all = np.arctan2(y_all, x_all)

        positive_x = x_all > 0.0
        if not np.any(positive_x):
            raise RuntimeError("No element centroids with x > 0 were found")

        positive_indices = np.where(positive_x)[0]
        target_index = positive_indices[np.argmin(np.abs(theta_all[positive_indices]))]
        theta_target = theta_all[target_index]

        angle_difference = np.arctan2(
            np.sin(theta_all - theta_target),
            np.cos(theta_all - theta_target),
        )

        row_mask = positive_x & (np.abs(angle_difference) < 1.0e-10)

        if np.count_nonzero(row_mask) < 2:
            row_mask = positive_x & np.isclose(
                theta_all,
                theta_target,
                rtol=0.0,
                atol=1.0e-8,
            )

        profile = current.loc[
            row_mask,
            ["x", "y", "vel_abs_x", "vel_abs_y", "pressure"],
        ].copy()

        if len(profile) < 2:
            raise RuntimeError(
                "Could not identify a complete radial row of element centroids"
            )

        profile["r"] = np.hypot(profile["x"], profile["y"])
        profile["theta_deg"] = np.degrees(np.arctan2(profile["y"], profile["x"]))
        profile = profile.sort_values("r").reset_index(drop=True)

        self.radius = profile["r"].to_numpy(dtype=float)

        # ------------------------------------------------------------------
        # Numerical values.
        #
        # The selected row is essentially the +x radial line, so the Cartesian
        # reconstructed absolute-frame components are used directly:
        #
        #   vel_abs_y -> u_theta
        #   vel_abs_x -> u_r
        # ------------------------------------------------------------------
        self.u_theta_numerical = profile["vel_abs_y"].to_numpy(dtype=float)
        self.u_r_numerical = profile["vel_abs_x"].to_numpy(dtype=float)
        self.pressure_numerical = profile["pressure"].to_numpy(dtype=float)

        # ------------------------------------------------------------------
        # Analytical Taylor-Couette solution evaluated at the exact same
        # centroid radii.
        # ------------------------------------------------------------------
        self.u_theta_analytical = self._velocity_solution(self.radius)
        self.u_r_analytical = np.zeros_like(self.radius)

        pressure_analytical_raw = self._pressure_solution(self.radius)

        # Pressure is defined only up to an additive constant. Shift numerical
        # and analytical pressure using the sampled centroid closest to the
        # pressure-pin radius.
        pressure_reference_index = np.argmin(
            np.abs(self.radius - self.PRESSURE_REFERENCE_RADIUS)
        )

        self.pressure_numerical_shifted = (
            self.pressure_numerical - self.pressure_numerical[pressure_reference_index]
        )

        self.pressure_analytical_shifted = (
            pressure_analytical_raw - pressure_analytical_raw[pressure_reference_index]
        )

        # ------------------------------------------------------------------
        # Validation metric for tangential velocity.
        #
        # Normalize by the physical inner-wall speed |Omega_i R_i|.
        # ------------------------------------------------------------------
        wall_speed = abs(self.OMEGA_I * self.RI)

        if wall_speed == 0.0:
            raise RuntimeError("Taylor-Couette wall-speed normalization is zero")

        self.u_theta_error = np.abs(self.u_theta_numerical - self.u_theta_analytical)

        self.u_theta_error_normalized = self.u_theta_error / wall_speed

        # Additional quantities retained for inspection.
        self.u_r_normalized = np.abs(self.u_r_numerical) / wall_speed

        self.pressure_error = np.abs(
            self.pressure_numerical_shifted - self.pressure_analytical_shifted
        )

        pressure_span = np.max(self.pressure_analytical_shifted) - np.min(
            self.pressure_analytical_shifted
        )

        if pressure_span == 0.0:
            raise RuntimeError("Taylor-Couette analytical pressure span is zero")

        self.pressure_error_normalized = self.pressure_error / pressure_span

        lower_bound = float(self.getParam("validation_lower_bound"))
        upper_bound = float(self.getParam("validation_upper_bound"))
        pressure_upper_bound = float(self.getParam("pressure_validation_upper_bound"))

        self.min_error = lower_bound * np.ones_like(self.u_theta_error_normalized)

        self.max_error = upper_bound * np.ones_like(self.u_theta_error_normalized)

        self.pressure_min_error = np.zeros_like(self.pressure_error_normalized)

        self.pressure_max_error = pressure_upper_bound * np.ones_like(
            self.pressure_error_normalized
        )

        # ------------------------------------------------------------------
        # Write a compact comparison CSV for review.
        # ------------------------------------------------------------------
        comparison = pd.DataFrame(
            {
                "r": self.radius,
                "x": profile["x"].to_numpy(dtype=float),
                "y": profile["y"].to_numpy(dtype=float),
                "theta_deg": profile["theta_deg"].to_numpy(dtype=float),
                "u_r_numerical": self.u_r_numerical,
                "u_r_analytical": self.u_r_analytical,
                "u_r_normalized": self.u_r_normalized,
                "u_theta_numerical": self.u_theta_numerical,
                "u_theta_analytical": self.u_theta_analytical,
                "u_theta_error": self.u_theta_error,
                "u_theta_error_normalized": self.u_theta_error_normalized,
                "pressure_numerical": self.pressure_numerical,
                "pressure_numerical_shifted": self.pressure_numerical_shifted,
                "pressure_analytical_shifted": self.pressure_analytical_shifted,
                "pressure_error": self.pressure_error,
                "pressure_error_normalized": self.pressure_error_normalized,
            }
        )

        comparison.to_csv(
            "taylor_couette_analytical_comparison.csv",
            index=False,
        )

        print(
            "Taylor-Couette centroid-row angle: " f"{np.degrees(theta_target):.12e} deg"
        )

        print(
            "Maximum normalized tangential-velocity error: "
            f"{np.max(self.u_theta_error_normalized):.12e}"
        )

        print(
            "RMS normalized tangential-velocity error: "
            f"{np.sqrt(np.mean(self.u_theta_error_normalized**2)):.12e}"
        )

        print(
            "Maximum normalized pressure-profile error: "
            f"{np.max(self.pressure_error_normalized):.12e}"
        )

        print(
            "RMS normalized pressure-profile error: "
            f"{np.sqrt(np.mean(self.pressure_error_normalized**2)):.12e}"
        )

    @classmethod
    def _velocity_coefficients(cls):
        """Return A and B for u_theta(r) = A r + B/r."""
        denominator = cls.RO**2 - cls.RI**2

        coefficient_a = (
            cls.OMEGA_O * cls.RO**2 - cls.OMEGA_I * cls.RI**2
        ) / denominator

        coefficient_b = (
            cls.RI**2 * cls.RO**2 * (cls.OMEGA_I - cls.OMEGA_O) / denominator
        )

        return coefficient_a, coefficient_b

    @classmethod
    def _velocity_solution(cls, radius):
        """Steady laminar Taylor-Couette tangential velocity."""
        coefficient_a, coefficient_b = cls._velocity_coefficients()

        return coefficient_a * radius + coefficient_b / radius

    @classmethod
    def _pressure_solution(cls, radius):
        """
        Analytical pressure distribution apart from an arbitrary constant.

        Radial equilibrium gives

            dp/dr = rho * u_theta^2 / r,

        with

            u_theta = A r + B/r.

        Therefore,

            p(r) = rho [
                (A^2/2) r^2
                + 2 A B ln(r)
                - B^2/(2 r^2)
            ] + C.
        """
        coefficient_a, coefficient_b = cls._velocity_coefficients()

        return cls.RHO * (
            0.5 * coefficient_a**2 * radius**2
            + 2.0 * coefficient_a * coefficient_b * np.log(radius)
            - 0.5 * coefficient_b**2 / radius**2
        )

    @staticmethod
    def _find_centroid_csv():
        """Locate the ElementValueSampler CSV produced by the relative case."""
        patterns = [
            "taylor_couette_2d_rel_out_centroid_profile_*.csv",
            "*centroid_profile_*.csv",
            "*centroid_profile*.csv",
        ]

        for pattern in patterns:
            candidates = sorted(Path(".").glob(pattern))
            if candidates:
                return candidates[-1]

        raise RuntimeError("Could not find the Taylor-Couette centroid-profile CSV")

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()

        params.addRequiredParam(
            "validation_lower_bound",
            "Lower bound for normalized tangential-velocity error",
        )

        params.addRequiredParam(
            "validation_upper_bound",
            "Upper bound for normalized tangential-velocity error",
        )

        params.addRequiredParam(
            "pressure_validation_upper_bound",
            "Upper bound for normalized pressure-profile error",
        )

        return params

    def testValidation(self):
        self.addVectorData(
            "u_theta_error",
            (self.radius, "Radius", "m"),
            (
                self.u_theta_error_normalized,
                "Normalized tangential-velocity error",
                "-",
            ),
            bounds=(self.min_error, self.max_error),
        )

        self.addVectorData(
            "pressure_error",
            (self.radius, "Radius", "m"),
            (
                self.pressure_error_normalized,
                "Normalized pressure-profile error",
                "-",
            ),
            bounds=(
                self.pressure_min_error,
                self.pressure_max_error,
            ),
        )
