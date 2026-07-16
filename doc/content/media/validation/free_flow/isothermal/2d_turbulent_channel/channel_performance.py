#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from TestHarness.resultsstore.reader import ResultsReader
from TestHarness.resultsstore.testdatafilters import TestDataFilter
from TestHarness.resultsstore.utils import TestName

os.chdir(os.path.dirname(os.path.realpath(__file__)))

test_names = {
    r"$Re_\tau = 395$": "channel_395",
    r"$Re_\tau = 590$": "channel_590",
}

# Load the timing history for both turbulent channel validation cases.
with ResultsReader("civet_tests_open_pronghorn_validation") as ctx:
    collection = ctx.reader.get_latest_push_results(50)
    test_results = {
        label: collection.get_tests(
            TestName("free_flow/isothermal/ercoftac_032_channel_flow", test_name),
            (TestDataFilter.STATUS, TestDataFilter.TIMING),
        )
        for label, test_name in test_names.items()
    }


def timing_dataframe(tests):
    """Return chronologically ordered timing data for non-timeout runs."""
    tests = [test for test in tests if test.status_value != "TIMEOUT"]
    df = pd.DataFrame(
        {
            "Date": [test.result.time for test in tests],
            "Runtime": [test.run_time for test in tests],
        }
    )
    df["Date"] = pd.to_datetime(df["Date"], errors="coerce")
    df = df.dropna(subset=["Date", "Runtime"]).sort_values(by="Date")
    df["5pt Avg"] = df["Runtime"].rolling(window=5, min_periods=1).mean()
    return df


timings = {label: timing_dataframe(tests) for label, tests in test_results.items()}
dates = sorted(
    set().union(*(set(df["Date"]) for df in timings.values() if not df.empty))
)
date_positions = {date: position for position, date in enumerate(dates)}

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

if not dates:
    ax.text(
        0.5,
        0.5,
        "No valid timing data available to plot.",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=12,
        transform=ax.transAxes,
    )
else:
    colors = ("#1f1f1f", "#0072B2")
    for (label, df), color in zip(timings.items(), colors):
        if df.empty:
            continue

        x_pos = np.array([date_positions[date] for date in df["Date"]])
        ax.plot(
            x_pos,
            df["Runtime"],
            marker="o",
            color=color,
            label=f"{label} Runtime",
        )
        ax.plot(
            x_pos,
            df["5pt Avg"],
            linestyle="--",
            color=color,
            label=f"{label} 5-point Avg",
        )

    ax.set_xticks(np.arange(len(dates)))
    ax.set_xticklabels([date.strftime("%Y-%m-%d") for date in dates], rotation=90)
    ax.set_xlabel("Date")
    ax.set_ylabel("Runtime (seconds)")
    ax.legend()
    ax.grid(True, color="gray", linestyle="--", linewidth=0.5, alpha=0.4)

fig.subplots_adjust(top=0.9, bottom=0.25, left=0.10)

plt.savefig("channel_performance.png", dpi=300)
plt.close()
