import unittest
from types import SimpleNamespace

import numpy as np

from pyPPG.fiducials import FpCollection

if not hasattr(np, "NaN"):
    np.NaN = np.nan


class PeakOnsetWindowBoundaryTests(unittest.TestCase):
    @staticmethod
    def detect(samples):
        sampling_rate = 100
        time = np.arange(samples, dtype=float) / sampling_rate
        ppg = np.sin(2 * np.pi * 1.2 * time)
        detector = FpCollection(SimpleNamespace(ppg=ppg, fs=sampling_rate))
        return detector.get_peak_onset()

    def test_accepts_exact_ten_second_window_at_100_hz(self):
        peaks, onsets = self.detect(1000)

        self.assertGreater(len(peaks), 0)
        self.assertGreater(len(onsets), 0)

    def test_rejects_input_shorter_than_detector_window(self):
        with self.assertRaisesRegex(ValueError, "requires at least 750 samples"):
            self.detect(999)


if __name__ == "__main__":
    unittest.main()
