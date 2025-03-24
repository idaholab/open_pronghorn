from strouhal import compute_strouhal
import unittest
import mooseutils
import os

class TestCylinder(unittest.TestCase):
    def test(self):
        executable = mooseutils.find_moose_executable_recursive(os.getcwd())
        cli_args = ['-i'] + ['flow.i']
        out = mooseutils.run_executable(executable, *cli_args, mpi=24, suppress_output=False)
        St = compute_strouhal("flow_csv.csv")
        self.assertTrue(St < 0.3050 and St > 0.2950)
