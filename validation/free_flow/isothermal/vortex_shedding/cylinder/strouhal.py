#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os

# Function responsible for doing the validation cehck on the Strouhal number
def run_validation_check(filenames):

  os.chdir(os.path.dirname(os.path.realpath(__file__)))
  d_cylinder = 0.1
  u_bulk = 1.0

  df  = []
  signals = []
  for filename in filenames:
    df.append(pd.read_csv(filename))
    signals.append(df['lift_coeff'])

  t = df[0]['time']
  St = []

  for signal in signals:

    fft_result = np.fft.fft(signal)

    frequencies = np.fft.fftfreq(len(fft_result), d=0.001)
    magnitude = np.abs(fft_result)

    peak_index = np.argmax(magnitude[:len(magnitude)//2])
    vortex_shedding_frequency = frequencies[peak_index]

    St.append(d_cylinder * vortex_shedding_frequency / u_bulk)

  st_df = pd.DataFrame(St, columns=['Strouhal'])
  st_df.to_csv('strouhal.csv', index=False)

  st_in_bounds = True
  for number in St:
    st_in_bounds = (number > 0.295 and number < 0.305)
    if (st_in_bounds):
      break

  return st_in_bounds
