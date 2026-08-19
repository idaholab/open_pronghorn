#!/usr/bin/env python3
import os
import subprocess
import sys

import numpy as np
import pandas as pd


def csv_matches_gold(current, gold, rel_err, abs_zero):
    if list(current.columns) != list(gold.columns) or len(current) != len(gold):
        return False

    for column in current.columns:
        current_values = current[column]
        gold_values = gold[column]
        current_numeric = pd.to_numeric(current_values, errors="coerce")
        gold_numeric = pd.to_numeric(gold_values, errors="coerce")

        if current_numeric.notna().all() and gold_numeric.notna().all():
            current_array = current_numeric.to_numpy(dtype=float, copy=True)
            gold_array = gold_numeric.to_numpy(dtype=float, copy=True)
            current_array[np.abs(current_array) < abs_zero] = 0.0
            gold_array[np.abs(gold_array) < abs_zero] = 0.0
            if not np.allclose(
                current_array,
                gold_array,
                rtol=rel_err,
                atol=abs_zero,
                equal_nan=False,
            ):
                return False
        elif not current_values.equals(gold_values):
            return False

    return True


def custom_evaluation(output):
    test_dir = os.path.dirname(os.path.realpath(__file__))
    command = [sys.executable, os.path.join(test_dir, "compute_halfwidths.py")]
    try:
        subprocess.run(command, cwd=test_dir, check=True)
    except subprocess.CalledProcessError:
        return False

    current = pd.read_csv(os.path.join(test_dir, "combined_halfwidths.csv"))
    gold = pd.read_csv(os.path.join(test_dir, "gold", "combined_halfwidths.csv"))
    if not csv_matches_gold(current, gold, rel_err=1e-3, abs_zero=1e-10):
        print("combined_halfwidths.csv does not match gold/combined_halfwidths.csv")
        return False

    return True
