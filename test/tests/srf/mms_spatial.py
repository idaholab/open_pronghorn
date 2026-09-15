import mms
import unittest
from mooseutils import fuzzyEqual, fuzzyAbsoluteEqual


class TestRotatingLid(unittest.TestCase):
    def test(self):
        velocity_labels = ["L2u", "L2v"]
        pressure_labels = ["L2p"]
        labels = velocity_labels + pressure_labels
        df1 = mms.run_spatial(
            "rotating_lid.i",
            6,
            y_pp=labels,
            file_base="rotating_lid",
        )

        fig = mms.ConvergencePlot(xlabel="Element Size ($h$)", ylabel="$L_2$ Error")
        fig.plot(
            df1,
            label=labels,
            marker="o",
            markersize=8,
            num_fitted_points=3,
            slope_precision=1,
        )
        fig.save("rotating_cavity.png")
        for key, value in fig.label_to_slope.items():
            print("%s, %f" % (key, value))
            if key in velocity_labels:
                self.assertTrue(fuzzyAbsoluteEqual(value, 3.5, 0.25))
            else:
                self.assertTrue(fuzzyAbsoluteEqual(value, 1.0, 0.25))
