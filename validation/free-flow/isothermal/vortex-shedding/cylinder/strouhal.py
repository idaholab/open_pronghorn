#!/usr/bin/env python3

import numpy as np
import pandas as pd
import mooseutils
import os

# Function responsible for doing the validation cehck on the Strouhal number
def run_validation_check():

  df = pd.read_csv("flow_csv.csv")
  t = df['time']
  signal = df['lift_coeff']

  d_cylinder = 0.1
  u_bulk = 1.0

  fft_result = np.fft.fft(signal)

  frequencies = np.fft.fftfreq(len(fft_result), d=0.001)
  magnitude = np.abs(fft_result)

  peak_index = np.argmax(magnitude[:len(magnitude)//2])
  vortex_shedding_frequency = frequencies[peak_index]

  St = d_cylinder * vortex_shedding_frequency / u_bulk

  df = pd.DataFrame([St], columns=['Strouhal'])
  df.to_csv('strouhal.csv', index=False)

  st_in_bounds = (St > 0.295 and St < 0.305)

  if not st_in_bounds:
    raise RuntimeError("St number out of bounds!")

  return 0
