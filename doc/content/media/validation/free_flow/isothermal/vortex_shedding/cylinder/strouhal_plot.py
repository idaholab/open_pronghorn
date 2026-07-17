#!/usr/bin/env python3

import os

import pandas as pd

from benchmark_plot import plot_benchmark_quantity

SCRIPT_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
os.chdir(SCRIPT_DIRECTORY)

csv_file = "../../../../../../../../validation/free_flow/isothermal/vortex_shedding/cylinder/gold/strouhal.csv"
value = pd.read_csv(csv_file).iloc[0, 0]

plot_benchmark_quantity(
    value=value,
    bounds=(0.295, 0.305),
    title="Strouhal Number",
    symbol="St",
    output_file="strouhal.png",
)
