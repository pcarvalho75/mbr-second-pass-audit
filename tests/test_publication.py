"""Unit tests of the published-workbook layer: display digits, typed cells and the frozen reading rule."""
import math
import unittest

from _paths import REPO  # noqa: F401
from mbraudit import workbook as W


class Display(unittest.TestCase):
    def test_rules_write_the_published_digits(self):
        self.assertEqual(W.display(0.029730617608409986, 'g4'), '0.02973')
        self.assertEqual(W.display('0.029730617608409986 [0.01600110488867993, 0.030314691996332688]', 'g4'), '0.02973 [0.016, 0.03031]')
        self.assertEqual(W.display('0.5374129431034861 to 2.239738690751491', 'f3'), '0.537 to 2.240')
        self.assertEqual(W.display('unavailable', 'f3'), 'unavailable')
        self.assertEqual(W.display(17023, 'thousands'), '17,023')
        self.assertEqual(W.display('[25.933417857691722%, 34.043655086749204%]', 'percent_or_f4'), '[25.93%, 34.04%]')
        self.assertEqual(W.display('[0.021608503447771987, 0.06389229461481327]', 'percent_or_f4'), '[0.0216, 0.0639]')
        self.assertEqual(W.display(0.20422732841879032), '0.20422732841879032')

    def test_cells_display_their_text(self):
        for text in ['292194', '17,023', '0.0666', '0.00066612239592853', '0.20422732841879032', '-0.0620', 'H6-1', '-0.0', '05']:
            value, fmt = W.to_cell(text)
            self.assertEqual(W.shown(value, fmt) if fmt else value, text)
        self.assertEqual(W.to_cell('0.20422732841879032'), ('0.20422732841879032', None))  # 17 digits stay text
        self.assertEqual(W.to_cell('17,023'), (17023, '#,##0'))


class Reading(unittest.TestCase):
    def test_frozen_reading_rule(self):
        row = {'endpoint': 'T1', 'numerator': 'T', 'fraction_ratio': 1.153, 'low95': 0.537, 'high95': 2.24}
        self.assertEqual(W.reading(row), 'Higher observed wrong-organism fraction in second-pass-only events. Interval includes no difference.')
        row = {'endpoint': 'T1', 'numerator': 'L', 'fraction_ratio': 3.027, 'low95': 1.561, 'high95': 6.196}
        self.assertEqual(W.reading(row), 'Above one: aligns directionally with Astral L/R. Finite interval excludes one.')
        row = {'endpoint': 'T2', 'numerator': 'positive_60_100', 'fraction_ratio': 0.239, 'low95': math.nan, 'high95': math.nan}
        self.assertTrue(W.reading(row).endswith('Interval unavailable or infinite; finite interval interpretation unavailable.'))


if __name__ == '__main__':
    unittest.main()
