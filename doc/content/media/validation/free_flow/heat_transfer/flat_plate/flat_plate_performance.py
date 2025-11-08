#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from TestHarness.resultsstore.reader import ResultsReader
import os

os.chdir(os.path.dirname(os.path.realpath(__file__)))

# We load the results
reader = ResultsReader("civet_tests_open_pronghorn_validation")
results = reader.getTestResults("free_flow/heat_transfer/flat_plate","flat_plate")

df = pd.DataFrame({
    'Date': [v.results.time for v in results],
    'Runtime': [v.run_time for v in results]
})

df = df.sort_values(by='Date')
df['DateStr'] = df['Date'].dt.strftime('%Y-%m-%d')
df['5pt Avg'] = df['Runtime'].rolling(window=5, min_periods=1).mean()
x_pos = np.arange(len(df))

# Formal appearance
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.titlesize": 14,
    "axes.titleweight": 'bold'
})

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x_pos, df['Runtime'], marker='o', color='#1f1f1f', label='Runtime')                # black
ax.plot(x_pos, df['5pt Avg'], marker='^', linestyle='-.', color='#D55E00', label='5-point Avg')  # reddish-orange

# X-axis
ax.set_xticks(x_pos)
ax.set_xticklabels(df['DateStr'], rotation=90)

# Axis labels only (no title)
ax.set_xlabel("Date")
ax.set_ylabel("Runtime (seconds)")
ax.legend()
ax.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.4)

fig.subplots_adjust(top=0.9, bottom=0.25, left=0.10)

# Export to PNG
output_file = "tlat_plate_performance.png"
plt.savefig(output_file, dpi=300)
plt.close()
