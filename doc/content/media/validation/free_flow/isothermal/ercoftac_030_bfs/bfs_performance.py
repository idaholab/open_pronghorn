#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from TestHarness.resultsstore.reader import ResultsReader
from TestHarness.resultsstore.testdatafilters import TestDataFilter
from TestHarness.resultsstore.utils import TestName
import os

# To make the relative paths work
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Load results
with ResultsReader("civet_tests_open_pronghorn_validation") as ctx:
    collection = ctx.reader.get_latest_push_results(50)
    tests = collection.get_tests(
        TestName("free_flow/isothermal/ercoftac_030_bfs", "bfs_030"),
        (TestDataFilter.STATUS, TestDataFilter.TIMING),
    )

# Get data for tests that didn't timeout
tests = [v for v in tests if v.status_value != "TIMEOUT"]
df = pd.DataFrame(
    {"Date": [v.result.time for v in tests], "Runtime": [v.run_time for v in tests]}
)

# Convert 'Date' to datetime, coercing errors to NaT
df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
df = df.sort_values(by="Date")

# Drop rows with NaT in 'Date'
df = df.dropna(subset=["Date"])

# Formal appearance
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 12,
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
    }
)

fig, ax = plt.subplots(figsize=(10, 6))

if df.empty:
    # Plot an empty figure with a message
    ax.text(
        0.5,
        0.5,
        "No valid datetime data available to plot.",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=12,
        transform=ax.transAxes,
    )
else:
    df["DateStr"] = df["Date"].dt.strftime("%Y-%m-%d")
    df["5pt Avg"] = df["Runtime"].rolling(window=5, min_periods=1).mean()
    x_pos = np.arange(len(df))

    # Plot data
    ax.plot(x_pos, df["Runtime"], marker="o", color="#1f1f1f", label="Runtime")  # black
    ax.plot(
        x_pos,
        df["5pt Avg"],
        marker="^",
        linestyle="-.",
        color="#D55E00",
        label="5-point Avg",
    )  # reddish-orange

    # X-axis
    ax.set_xticks(x_pos)
    ax.set_xticklabels(df["DateStr"], rotation=90)

    # Axis labels only (no title)
    ax.set_xlabel("Date")
    ax.set_ylabel("Runtime (seconds)")
    ax.legend()
    ax.grid(True, color="gray", linestyle="--", linewidth=0.5, alpha=0.4)

fig.subplots_adjust(top=0.9, bottom=0.25, left=0.10)

# Export to PNG
output_file = "bfs_performance.png"
plt.savefig(output_file, dpi=300)
plt.show()
