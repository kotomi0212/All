import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np
from balance_analysis import read_recording, analyze, mode_metrics


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name)/'test.csv'

    def fixture(self, times=None, nan=False, duplicate_header=True):
        t = np.arange(4000)/100 if times is None else times
        a = np.sin(2*np.pi*t)
        h = 2*np.cos(2*np.pi*t)
        with self.path.open('w', encoding='cp932', newline='') as f:
            w = csv.writer(f)
            header = ['', 'xix1', 'xix2', 'xiy1', 'xiy2', 'axhd', 'ayhd']
            w.writerow(['記録周波数(Hz)', '100'])
            if duplicate_header:
                w.writerow(header)
                w.writerow(['単位']+['mm/s/s']*6)
                w.writerow(['MIN']+[0]*6)
            w.writerow(header)
            for i in range(len(t)):
                w.writerow([t[i], 'nan' if nan and i == 10 else a[i], h[i], a[i], h[i], a[i]+h[i], a[i]+h[i]])
        return self.path

    def test_pure_ankle_and_zero(self):
        a = np.array([-1, 0, 1])
        self.assertAlmostEqual(mode_metrics(a, np.zeros(3))['ankle_percent'], 100)
        self.assertAlmostEqual(mode_metrics(np.zeros(3), a)['hip_percent'], 100)
        self.assertIsNone(mode_metrics(np.zeros(3), a)['rxi'])
        self.assertIsNone(mode_metrics(np.ones(3), np.ones(3))['ankle_percent'])

    def test_rms_ratio_not_variance_share_and_offset_invariance(self):
        a = np.array([-1, 0, 1])
        m = mode_metrics(a + 100, -2*a + 200)
        self.assertAlmostEqual(m['ankle_percent'], 100/3)
        self.assertAlmostEqual(m['hip_percent'], 200/3)
        self.assertAlmostEqual(m['rxi'], 2)
        self.assertAlmostEqual(m['cxi'], -1)

    def test_repeated_header_and_interval(self):
        r = read_recording(self.fixture())
        out = analyze(r)
        self.assertEqual(out['samples'], 3000)
        self.assertAlmostEqual(out['axes']['AP']['rxi'], 2)
        self.assertLess(abs(out['axes']['AP']['cxi']), 1e-12)
        self.assertEqual(out['warnings'], [])

    def test_reconstruction_warning(self):
        r = read_recording(self.fixture())
        r.signals['ayhd'] *= 2
        self.assertTrue(any(w.startswith('ML:') for w in analyze(r)['warnings']))

    def test_bad_values_and_missing_samples(self):
        with self.assertRaisesRegex(ValueError, '欠損'):
            read_recording(self.fixture(nan=True))
        with self.assertRaisesRegex(ValueError, '不均一'):
            read_recording(self.fixture(times=np.delete(np.arange(4000)/100, 10)))

    def test_duplicate_time_and_bad_interval(self):
        t = np.arange(4000)/100
        t[10] = t[9]
        with self.assertRaisesRegex(ValueError, '昇順'):
            read_recording(self.fixture(times=t))
        r = read_recording(self.fixture())
        for start, end in [(35, 5), (-1, 35), (5, 41), (float('nan'), 35)]:
            with self.assertRaises(ValueError):
                analyze(r, start, end)

    def test_raw_only_rejected(self):
        self.path.write_text('time,Fx,Fy,AccX1\n0,1,1,1', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, '処理済みモード'):
            read_recording(self.path)


if __name__ == '__main__':
    unittest.main()
