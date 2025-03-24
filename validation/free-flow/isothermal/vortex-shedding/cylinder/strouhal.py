import numpy as np
import pandas as pd

def compute_strouhal(filename):

  df = pd.read_csv(filename)
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

  print(St)

  return St
