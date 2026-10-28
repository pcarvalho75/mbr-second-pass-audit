"""The post hoc carryover check on a hand-counted synthetic search: two human-only runs, one yeast-only and one
E. coli-only source run acquired before them, and one blank."""
import unittest

import numpy as np
import pandas as pd

from _paths import REPO  # noqa: F401
from mbraudit.pg250.carryover import analyse

DESIGN = pd.DataFrame({
    'sample': ['SYN_E1', 'SYN_Y1', 'SYN_B1', 'SYN_H1', 'SYN_H2'],
    'kind': ['single', 'single', 'blank', 'single', 'single'],
    'present_species': ['E', 'Y', '', 'H', 'H'],
    'filename_date': ['20250312'] * 5,
    'acquisition_id': [1, 3, 4, 5, 6]})


def frame(stage, rows):
    out = pd.DataFrame(rows, columns=['sample', 'Modified.Sequence', 'Precursor.Charge', 'species', 'Precursor.Quantity'])
    out['mask'], out['stage'] = 'F1', stage
    return out


FIRST = frame('first', [('SYN_H1', 'HHH', 2, 'H', 50.), ('SYN_H1', 'YYYA', 2, 'Y', 20.), ('SYN_H2', 'HHH', 2, 'H', 60.),
                        ('SYN_Y1', 'YYYA', 2, 'Y', 100.), ('SYN_Y1', 'YYYB', 2, 'Y', 1000.), ('SYN_E1', 'EEEA', 2, 'E', 50.)])
FINAL = frame('final', [('SYN_H1', 'HHH', 2, 'H', 50.), ('SYN_H1', 'YYYA', 2, 'Y', 20.), ('SYN_H1', 'YYYB', 2, 'Y', 100.),
                        ('SYN_H1', 'EEEA', 2, 'E', 5.), ('SYN_H1', 'HHHB', 2, 'H', 7.),
                        ('SYN_H2', 'HHH', 2, 'H', 60.), ('SYN_H2', 'YYYB', 2, 'Y', 10.), ('SYN_H2', 'HHHC', 2, 'H', 8.),
                        ('SYN_Y1', 'YYYA', 2, 'Y', 100.), ('SYN_Y1', 'YYYB', 2, 'Y', 1000.), ('SYN_E1', 'EEEA', 2, 'E', 50.),
                        ('SYN_B1', 'HHH', 2, 'H', 3.), ('SYN_B1', 'YYYB', 2, 'Y', 4.)])


class Carryover(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        t = analyse(FIRST, FINAL, DESIGN)
        cls.value = {(r.section, r.sample, r.organism, r.measure): r.value for r in t.itertuples()}

    def get(self, section, measure, sample='', organism=''):
        return self.value[(section, sample, organism, measure)]

    def test_classes_and_run_order(self):
        self.assertEqual((self.get('check', 'added wrong-organism'), self.get('check', 'added species-exclusive')), (3, 5))
        self.assertEqual((self.get('check', 'retained wrong-organism'), self.get('check', 'retained species-exclusive')), (1, 3))
        self.assertEqual(self.get('run order', 'added yeast', 'SYN_H1', 'Y'), 1)
        self.assertEqual(self.get('run order', 'added E. coli', 'SYN_H1', 'E'), 1)
        self.assertAlmostEqual(self.get('run order', 'added wrong-organism fraction, %', 'SYN_H1'), 200 / 3)
        self.assertAlmostEqual(self.get('run order', 'added wrong-organism fraction, %', 'SYN_H2'), 50.0)
        self.assertEqual((self.get('organism', 'added yeast', organism='Y'), self.get('organism', 'added E. coli', organism='E')), (2, 1))

    def test_intensity_against_the_source_runs(self):
        # yeast: 100/1000 in SYN_H1 and 10/1000 in SYN_H2, log10 ratios -1 and -2
        self.assertAlmostEqual(self.get('intensity', 'median ratio to the source median', organism='Y'), 10 ** -1.5)
        self.assertAlmostEqual(self.get('intensity', 'interquartile range of log10 ratios', organism='Y'), 0.5)
        self.assertEqual(self.get('intensity', 'added yeast ions accepted in the source runs', organism='Y'), 2)
        self.assertAlmostEqual(self.get('intensity', 'median ratio to the source median', organism='E'), 0.1)

    def test_recurrence_and_blanks(self):
        self.assertEqual(self.get('recurrence', 'distinct wrong-organism precursors'), 2)
        self.assertEqual(self.get('recurrence', 'added in 1 of the six runs'), 1)
        self.assertEqual(self.get('recurrence', 'added in 2 of the six runs'), 1)
        self.assertEqual(self.get('blanks', 'first accepted target precursors', 'SYN_B1'), 0)
        self.assertEqual((self.get('blanks', 'final accepted H', 'SYN_B1', 'H'), self.get('blanks', 'final accepted Y', 'SYN_B1', 'Y')), (1, 1))


if __name__ == '__main__':
    unittest.main()
