#!/usr/bin/env python3
import os
import subprocess
import sys


def custom_evaluation(output):
    test_dir = os.path.dirname(os.path.realpath(__file__))
    command = [sys.executable, os.path.join(test_dir, "compute_halfwidths.py")]
    try:
        subprocess.run(command, cwd=test_dir, check=True)
    except subprocess.CalledProcessError:
        return False
    return True
