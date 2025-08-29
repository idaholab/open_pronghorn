#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from TestHarness.resultsreader.reader import TestHarnessResultsReader

reader = TestHarnessResultsReader("civet_tests_open_pronghorn_validation")
results = reader.getTestResults("free_flow/isothermal/ercoftac_030_bfs_copy","bfs_030")

# We populate the dataframe
dates = []
runtimes = []
for c in results:
  dates.append(c.results.time)
  runtimes.append(c.run_time)

df = pd.DataFrame({
    'Date': dates,
    'Runtime': runtimes
})

df = df.sort_values(by='Date')
df['DateStr'] = df['Date'].dt.strftime('%Y-%m-%d')
df['3pt Avg'] = df['Runtime'].rolling(window=3, min_periods=1).mean()
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
fig, ax = plt.subplots(figsize=(15, 6))
ax.plot(x_pos, df['Runtime'], marker='o', color='#1f1f1f', label='Runtime')                # black
ax.plot(x_pos, df['3pt Avg'], marker='s', linestyle='--', color="#629ABB", label='3-point Avg')  # blue
ax.plot(x_pos, df['5pt Avg'], marker='^', linestyle='-.', color='#D55E00', label='5-point Avg')  # reddish-orange

# X-axis
ax.set_xticks(x_pos)
ax.set_xticklabels(df['DateStr'], rotation=90)

# Axis labels only (no title)
ax.set_xlabel("Date")
ax.set_ylabel("Runtime (seconds)")
ax.legend()
ax.grid(True, color='gray', linestyle='--', linewidth=0.5, alpha=0.4)

fig.subplots_adjust(left=0.3, right=0.7, top=0.7, bottom=0.3)

# Export to PNG
output_file = "bfs_performance.png"
plt.savefig(output_file, dpi=300)
plt.show()
