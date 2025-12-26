from TestHarness.validation import ValidationCase

import numpy as np
import pandas as pd


class TestCase(ValidationCase):
    """Validation for Chen aerosol drift-flux channel benchmark.

    Compares (i) axial velocity profiles and (ii) aerosol concentration profiles
    extracted from LineValueSampler outputs (x = 0.2, 0.4, 0.6 m) against
    digitized experimental data.

    The test follows the same pattern as the ERCOFTAC channel-flow validations:
      1) compute the current error-to-experiment
      2) compute a *reference* (gold) OpenPronghorn error-to-experiment
      3) require the current error to remain inside a configurable envelope
         around the reference error.

    Expected files in the test directory:
      - chen_steady_csv_sampler_line_x_0d2_0002.csv
      - chen_steady_csv_sampler_line_x_0d4_0002.csv
      - chen_steady_csv_sampler_line_x_0d6_0002.csv

    Expected files in gold/:
      - data_u_x_0d2.csv, data_u_x_0d4.csv, data_u_x_0d6.csv
      - data_c_x_0d2.csv, data_c_x_0d4.csv, data_c_x_0d6.csv
      - chen_steady_csv_sampler_line_x_0d2_0002.csv (reference)
      - chen_steady_csv_sampler_line_x_0d4_0002.csv (reference)
      - chen_steady_csv_sampler_line_x_0d6_0002.csv (reference)

    Experimental CSV convention:
      - col 1: u/Uref (or C/Cref)
      - col 2: z/H

    Sampler CSV convention:
      - columns include 'vel_x', 'aerosol', 'z' (meters)
    """

    # Keep in sync with chen_steady.i
    H = 0.4
    U_REF = 0.225  # inlet velocity in chen_steady.i
    C_REF = 1.0

    def initialize(self):
        self.lower_bound = float(self.getParam("validation_lower_bound"))
        self.upper_bound = float(self.getParam("validation_upper_bound"))

        cases = [
            ("0d2", 0.2),
            ("0d4", 0.4),
            ("0d6", 0.6),
        ]

        self._vector_payloads = []

        for tag, _xloc in cases:
            # Current and reference sampler outputs
            cur = self._load_sampler(f"chen_steady_csv_sampler_line_x_{tag}_0002.csv")
            ref = self._load_sampler(f"gold/chen_steady_csv_sampler_line_x_{tag}_0002.csv")

            # Experimental digitized data
            exp_u = self._load_exp_two_col(f"gold/data_u_x_{tag}.csv")
            exp_c = self._load_exp_two_col(f"gold/data_c_x_{tag}.csv")

            # Non-dimensional wall-normal coordinate
            zH_cur = cur["z"] / self.H
            zH_ref = ref["z"] / self.H

            # Interpolate experiment onto current sampler positions
            u_exp_on_cur = np.interp(zH_cur, exp_u["zH"], exp_u["val"])
            c_exp_on_cur = np.interp(zH_cur, exp_c["zH"], exp_c["val"])

            # Normalize current profiles
            u_cur = cur["vel_x"] / self.U_REF
            c_cur = cur["aerosol"] / self.C_REF

            # Normalize reference profiles, then interpolate them onto current z/H
            u_ref = ref["vel_x"] / self.U_REF
            c_ref = ref["aerosol"] / self.C_REF
            u_ref_on_cur = np.interp(zH_cur, zH_ref, u_ref)
            c_ref_on_cur = np.interp(zH_cur, zH_ref, c_ref)

            # Current and reference absolute errors to experiment
            u_err_cur = np.abs(u_cur - u_exp_on_cur)
            u_err_ref = np.abs(u_ref_on_cur - u_exp_on_cur)

            c_err_cur = np.abs(c_cur - c_exp_on_cur)
            c_err_ref = np.abs(c_ref_on_cur - c_exp_on_cur)

            # Envelope around reference error
            u_min = (1.0 + self.lower_bound) * u_err_ref
            u_max = (1.0 + self.upper_bound) * u_err_ref
            c_min = (1.0 + self.lower_bound) * c_err_ref
            c_max = (1.0 + self.upper_bound) * c_err_ref

            # Payloads (errors)
            self._vector_payloads.append(
                dict(
                    name=f"Ux_diff_x_{tag}",
                    x=zH_cur,
                    y=u_err_cur,
                    ybounds=(u_min, u_max),
                    xlabel="z/H",
                    ylabel="|u/Uref - u_exp|",
                )
            )
            self._vector_payloads.append(
                dict(
                    name=f"C_diff_x_{tag}",
                    x=zH_cur,
                    y=c_err_cur,
                    ybounds=(c_min, c_max),
                    xlabel="z/H",
                    ylabel="|C/Cref - C_exp|",
                )
            )

            # Payloads (profiles)
            self._vector_payloads.append(
                dict(
                    name=f"Ux_x_{tag}",
                    x=zH_cur,
                    y=u_cur,
                    ybounds=None,
                    xlabel="z/H",
                    ylabel="u/Uref",
                )
            )
            self._vector_payloads.append(
                dict(
                    name=f"C_x_{tag}",
                    x=zH_cur,
                    y=c_cur,
                    ybounds=None,
                    xlabel="z/H",
                    ylabel="C/Cref",
                )
            )

    @staticmethod
    def validParams():
        params = ValidationCase.validParams()
        params.addRequiredParam(
            "validation_lower_bound",
            "Lower bound for relative envelope around reference error-to-experiment.",
        )
        params.addRequiredParam(
            "validation_upper_bound",
            "Upper bound for relative envelope around reference error-to-experiment.",
        )
        return params

    def testValidation(self):
        for payload in self._vector_payloads:
            bounds = payload["ybounds"]
            if bounds is None:
                self.addVectorData(
                    payload["name"],
                    (payload["x"], payload["xlabel"], "-"),
                    (payload["y"], payload["ylabel"], "-"),
                )
            else:
                self.addVectorData(
                    payload["name"],
                    (payload["x"], payload["xlabel"], "-"),
                    (payload["y"], payload["ylabel"], "-"),
                    bounds=bounds,
                )

    @staticmethod
    def _load_exp_two_col(path):
        df = pd.read_csv(path)
        val = df[df.columns[0]].to_numpy(dtype=float)
        zH = df[df.columns[1]].to_numpy(dtype=float)
        idx = np.argsort(zH)
        return {"val": val[idx], "zH": zH[idx]}

    @staticmethod
    def _load_sampler(path):
        df = pd.read_csv(path, comment="#")
        required = {"vel_x", "aerosol", "z"}
        missing = required.difference(df.columns)
        if missing:
            raise RuntimeError(f"Sampler CSV '{path}' missing columns: {sorted(missing)}")
        df = df.sort_values("z")
        return {
            "vel_x": df["vel_x"].to_numpy(dtype=float),
            "aerosol": df["aerosol"].to_numpy(dtype=float),
            "z": df["z"].to_numpy(dtype=float),
        }
