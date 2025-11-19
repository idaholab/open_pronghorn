import pandas as pd
import numpy as np
from scipy.special import erfcx


def read_moose(path, independent_var):
    df = pd.read_csv(path)
    t_cols = [c for c in df.columns if c.startswith("T")]
    if len(t_cols) != 1:
        raise ValueError(
            f"Expected exactly one temperature column starting with 'T' in {path}, got {t_cols}"
        )
    df = df.rename(columns={t_cols[0]: "T"})
    # Keep standard columns if present; otherwise keep what exists
    keep_cols = [c for c in ["T", independent_var] if c in df.columns]
    return df[keep_cols].copy()


def analytic_approximate_solution_dht_erfcx(x, y):
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
        T_w = analytic_approximate_solution_dht_erfcx(x, 1e-12)
        return T_w + (T_w - T_b) * y / b

    return T_b + (T_inf - T_b) * F


def analytic_approximate_solution_bl(x, y):
    T_b = 600
    T_inf = 1000
    lambda_s = 0.2876
    lambda_f = 0.06808
    b = 0.01
    rho_f = 0.3525
    u_inf = 12.0
    cp = 1142.6
    mu = 3.95e-5

    Re_x = rho_f * x * u_inf / mu
    Pr = mu * cp / lambda_f
    delta = 4.64 / np.sqrt(Re_x) * x
    delta_T = delta * np.pow(13.0 / 14.0, 1.0 / 3.0) * np.pow(1.0 / Pr, 1.0 / 3.0)
    z = 3.0 / 2.0 * lambda_f / lambda_s * b / delta_T
    T_w = T_b + (T_inf - T_b) * z / (1 + z)

    if y <= 0:
        T_w = analytic_approximate_solution_bl(x, 1e-12)
        return T_w + (T_w - T_b) * y / b

    if y < delta_T and y > 0:
        return (
            T_w
            + 3.0 / 2.0 * (T_inf - T_w) / delta_T * y
            - 1.0 / 2.0 * (T_inf - T_w) / np.pow(delta_T, 3.0) * np.pow(y, 3.0)
        )

    return None
