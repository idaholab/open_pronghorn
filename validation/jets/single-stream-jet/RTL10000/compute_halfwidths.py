#!/usr/bin/env python3
"""
Stable half-width detection:
- halfwidth = smallest y>0 where vel_z < 0.5*Vc and stays below for min_consecutive points
- applies light smoothing (Gaussian-like convolution) to vel_z(y) before detection
- writes combined_halfwidths.csv
- deletes radial csv's after processing
"""

import os
import re
import numpy as np
import pandas as pd

OUTFILE = "combined_halfwidths.csv"

# Parameters
SMOOTH_WINDOW = 7
SMOOTH_PAD_MODE = "edge"
MIN_CONSECUTIVE = 5
MIN_SPAN_Y = None


def find_files():
    return [
        f
        for f in os.listdir(".")
        if os.path.isfile(f) and re.search(r"(?i)^jet_RTL10000_out_radial.*\.csv$", f)
    ]


def extract_index(fname):
    m = re.search(r"(?i)^jet_RTL10000_out_radial[_-]?0*([0-9]{1,3})", fname)
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
    ys = g["y"].to_numpy()
    vz = g["vel_z"].to_numpy()

    if len(ys) == 0:
        return np.nan

    mask0 = np.isclose(ys, 0.0)
    if mask0.any():
        return float(vz[mask0].mean())

    idx = np.argsort(ys)
    ys = ys[idx]
    vz = vz[idx]

    # If no sign change around y=0, use point closest to zero
    if ys[0] > 0 or ys[-1] < 0:
        return float(vz[np.argmin(np.abs(ys))])

    i = np.where(ys < 0)[0][-1]
    return float(vz[i] + (vz[i + 1] - vz[i]) * (0.0 - ys[i]) / (ys[i + 1] - ys[i]))


def smooth_signal(v, window=7, pad_mode="edge"):
    if window <= 1:
        return v.copy()
    if window % 2 == 0:
        window += 1

    x = np.linspace(-(window // 2), window // 2, window)
    sigma = max(window / 4, 1e-12)
    kernel = np.exp(-0.5 * (x / sigma) ** 2)
    kernel /= kernel.sum()

    pad = window // 2
    vpad = np.pad(v, pad, mode=pad_mode)
    return np.convolve(vpad, kernel, mode="valid")


def interp_y(v0, y0, v1, y1, vt):
    if np.isclose(v0, v1):
        return (y0 + y1) / 2.0
    t = (vt - v0) / (v1 - v0)
    return y0 + t * (y1 - y0)


def halfwidth_stable(g, vc, smooth_window=7, min_consec=5, min_span_y=None):
    if np.isnan(vc):
        return np.nan

    target = 0.5 * vc

    gpos = g[g["y"] > 0].sort_values("y").reset_index(drop=True)
    if len(gpos) < 2:
        return np.nan

    y = gpos["y"].to_numpy()
    v = gpos["vel_z"].to_numpy()
    vs = smooth_signal(v, window=smooth_window, pad_mode=SMOOTH_PAD_MODE)

    n = len(vs)
    for i in range(n - 1):
        if vs[i] >= target and vs[i + 1] < target:
            end = min(n, i + 1 + min_consec)
            if np.all(vs[i + 1 : end] < target):
                y_cross = interp_y(vs[i], y[i], vs[i + 1], y[i + 1], target)
                if min_span_y is not None:
                    span_y = y[end - 1] - y_cross
                    if span_y < min_span_y:
                        continue
                return float(y_cross)

    return np.nan


def main():
    files = find_files()
    if not files:
        print("No matching files found.")
        return

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
            g = df[df["z"] == z][["y", "vel_z"]].copy()
            vc = centerline_velocity(g)
            yh = halfwidth_stable(
                g,
                vc,
                smooth_window=SMOOTH_WINDOW,
                min_consec=MIN_CONSECUTIVE,
                min_span_y=MIN_SPAN_Y,
            )
            rows.append(
                [int(idx), fname, float(z), float(yh) if not np.isnan(yh) else np.nan]
            )

        # # Delete only after the whole file has been processed
        # os.remove(fname)

    out = pd.DataFrame(rows, columns=["file_index", "file_name", "z", "y_half"])
    out.to_csv(OUTFILE, index=False)
    print(f"Wrote {OUTFILE} ({len(out)} rows)")


if __name__ == "__main__":
    main()
