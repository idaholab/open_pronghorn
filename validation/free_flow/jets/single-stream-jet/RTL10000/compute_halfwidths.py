#!/usr/bin/env python3
"""
Half-width detection:
- halfwidth = interpolated y > 0 where smoothed vel_z first crosses below 0.5*Vc
- crossing is accepted only when smoothed vel_z stays below for min_consecutive points
- writes combined_halfwidths.csv
"""

import os
import re
import numpy as np
import pandas as pd

OUTFILE = "combined_halfwidths.csv"

SMOOTH_WINDOW = 7  # smoothing kernel size in samples (odd int >=3)
SMOOTH_PAD_MODE = "edge"  # mode for np.pad
MIN_CONSECUTIVE = (
    5  # require at least this many consecutive samples below threshold after crossing
)
MIN_SPAN_Y = None  # optionally require that span in y after crossing >= value (e.g. 0.005). None -> ignore


# ---------------- helpers ----------------
def find_files():
    return [
        f
        for f in os.listdir(".")
        if os.path.isfile(f) and re.search(r"(?i)^jet_RTL10000_out_radial_.*\.csv$", f)
    ]


def extract_index(fname):
    m = re.search(r"(?i)jet_RTL10000_out_radial[_-]?0*([0-9]{1,3})", fname)
    return int(m.group(1)) if m else None


def read_df(fname):
    df = pd.read_csv(fname)
    for c in ("vel_z", "y", "z"):
        if c not in df.columns:
            raise ValueError(f"{fname} missing column '{c}'")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df["z"] = pd.to_numeric(df["z"], errors="coerce")
    df["vel_z"] = pd.to_numeric(df["vel_z"], errors="coerce")
    return df.dropna(subset=["y", "z", "vel_z"])


def centerline_velocity(g):
    ys = g["y"].values
    vz = g["vel_z"].values
    if np.isclose(ys, 0.0).any():
        return float(vz[np.isclose(ys, 0.0)].mean())
    idx = np.argsort(ys)
    ys, vz = ys[idx], vz[idx]
    if ys[0] > 0 or ys[-1] < 0:
        return float(vz[np.argmin(np.abs(ys))])
    i = np.where(ys < 0)[0][-1]
    return float(vz[i] + (vz[i + 1] - vz[i]) * (0.0 - ys[i]) / (ys[i + 1] - ys[i]))


def smooth_signal(y, v, window=7, pad_mode="edge"):
    """Smooth the radial profile before threshold detection."""
    if window <= 1:
        return v.copy()
    if window % 2 == 0:
        window += 1
    kernel = np.exp(
        -0.5 * (np.linspace(-(window // 2), (window // 2), window) / (window / 4)) ** 2
    )
    kernel = kernel / kernel.sum()
    pad = window // 2
    vpad = np.pad(v, pad, mode=pad_mode)
    return np.convolve(vpad, kernel, mode="valid")


def interp_y(v0, y0, v1, y1, vt):
    if np.isclose(v0, v1):
        return (y0 + y1) / 2.0
    t = (vt - v0) / (v1 - v0)
    return y0 + t * (y1 - y0)


def halfwidth_stable(g, vc, smooth_window=7, min_consec=5, min_span_y=None):
    """
    Returns the first y > 0 where smoothed vel_z crosses below 0.5*vc and remains below
    for at least min_consec samples (and optional min_span_y in physical units).
    """
    if np.isnan(vc):
        return np.nan
    target = 0.5 * vc
    # restrict to positive y and sort
    gpos = g[g["y"] > 0].sort_values("y").reset_index(drop=True)
    if len(gpos) < 2:
        return np.nan
    y = gpos["y"].values
    v = gpos["vel_z"].values

    # Use smoothing only to make threshold detection insensitive to local oscillations.
    # The output remains a single interpolated half-width value, not a smoothed profile.
    vs = smooth_signal(y, v, window=smooth_window, pad_mode=SMOOTH_PAD_MODE)

    n = len(vs)
    # iterate over consecutive pairs (i,i+1) where the smoothed profile crosses below target
    for i in range(n - 1):
        if vs[i] >= target and vs[i + 1] < target:
            # check that vs[i+1 : i+1+min_consec] are all < target
            end = i + 1 + min_consec
            if end > n:
                continue
            if np.all(vs[i + 1 : end] < target):
                y_cross = interp_y(vs[i], y[i], vs[i + 1], y[i + 1], target)
                # optionally require that span from y_cross to y[end-1] >= min_span_y
                if min_span_y is not None:
                    span_y = y[end - 1] - y_cross if end - 1 >= 0 else 0.0
                    if span_y < min_span_y:
                        continue
                return float(y_cross)
    return np.nan


# ---------------- main ----------------
def main():
    files = find_files()
    if not files:
        raise SystemExit("No jet_RTL10000 radial output CSV files found")

    indexed = [(extract_index(f), f) for f in files]
    indexed = [(i, f) for i, f in indexed if i is not None]
    indexed.sort(key=lambda x: x[0])

    rows = []
    for idx, fname in indexed:
        try:
            df = read_df(fname)
        except Exception as e:
            print(f"Skipping {fname}: {e}")
            continue
        for z in np.sort(df["z"].unique()):
            g = df[df["z"] == z][["y", "vel_z"]]
            vc = centerline_velocity(g)
            yh = halfwidth_stable(
                g,
                vc,
                smooth_window=SMOOTH_WINDOW,
                min_consec=MIN_CONSECUTIVE,
                min_span_y=MIN_SPAN_Y,
            )
            rows.append([float(z), float(yh) if not np.isnan(yh) else np.nan])

    if not rows:
        raise SystemExit("No jet half-width rows were computed")

    out = pd.DataFrame(rows, columns=["z", "y_half"])
    if out["y_half"].isna().any():
        raise SystemExit("One or more jet half-width rows could not be computed")

    out.to_csv(OUTFILE, index=False)
    print(f"Wrote {OUTFILE} ({len(out)} rows)")


if __name__ == "__main__":
    main()
