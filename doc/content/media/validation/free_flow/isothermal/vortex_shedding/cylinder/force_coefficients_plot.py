#!/usr/bin/env python3

import os

import pandas as pd

from benchmark_plot import plot_benchmark_quantity

SCRIPT_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
os.chdir(SCRIPT_DIRECTORY)

csv_file = "../../../../../../../../validation/free_flow/isothermal/vortex_shedding/cylinder/gold/force_coefficients.csv"
data = pd.read_csv(csv_file).iloc[0]

quantities = (
    {
        "column": "Maximum drag coefficient",
        "title": "Maximum Drag Coefficient",
        "symbol": r"C_{D,\max}",
        "bounds": (3.22, 3.24),
        "output_file": "drag_coefficient.png",
    },
    {
        "column": "Maximum lift coefficient",
        "title": "Maximum Lift Coefficient",
        "symbol": r"C_{L,\max}",
        "bounds": (0.99, 1.01),
        "output_file": "lift_coefficient.png",
    },
)

for quantity in quantities:
    plot_benchmark_quantity(
        value=data[quantity["column"]],
        bounds=quantity["bounds"],
        title=quantity["title"],
        symbol=quantity["symbol"],
        output_file=quantity["output_file"],
    )
