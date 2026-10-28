"""Unit tests of the single-cell definitions: filters, species calls, classes, support, intervals, ratios."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from _paths import REPO  # noqa: F401
from mbraudit.filters import filter_mask
from mbraudit.species import SpeciesRegistry, accept_rows, classify, combined_factor, run_name
from mbraudit.singlecell import bootstrap as bb
from mbraudit.singlecell.accounting import collapse_labels, fraction_ratio, quantity_terciles, species_record, support_bin
from mbraudit.singlecell.composition import rate_ratio_ci
from mbraudit.singlecell.intervals import derive_counts
from mbraudit.singlecell.reports import group_units, natural
from mbraudit.singlecell.support import original_bin, wilson


class Filters(unittest.TestCase):
    def test_f1_f2_f3_are_literal_conjunctions(self):
        frame = pd.DataFrame({'Q.Value': [0.01, 0.01, 0.01, 0.02], 'Global.Q.Value': [0.01, 0.02, 0.01, 0.0],
                              'Lib.Q.Value': [0.0, 0.0, 0.01, 0.0], 'PG.Q.Value': [0.01, 0.0, 0.02, 0.0],
                              'Global.PG.Q.Value': [0.0, 0.0, 0.0, 0.0], 'Lib.PG.Q.Value': [0.0, 0.0, 0.0, 0.0]})
        self.assertEqual(filter_mask(frame, 'F1').tolist(), [True, False, True, False])
        self.assertEqual(filter_mask(frame, 'F2').tolist(), [True, False, False, False])
        self.assertEqual(filter_mask(frame, 'F3').tolist(), [True, True, True, False])
        with self.assertRaises(KeyError):
            filter_mask(frame, 'F9')

    def test_missing_fields_do_not_pass_once_filled_with_two(self):
        frame = pd.DataFrame({c: [np.nan] for c in ['Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value', 'Lib.PG.Q.Value']})
        self.assertFalse(filter_mask(frame.fillna(2.), 'F1').iloc[0])


class Species(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        root = Path(self.temp.name)
        (root / 'human.fasta').write_text('>sp|H1|H1_HUMAN\nMPEPTLDEK\n>sp|H2|H2_HUMAN\nMSAMPLEK\n', encoding='utf-8')
        (root / 'both.fasta').write_text('>sp|H1|H1_HUMAN\nMPEPTLDEK\n>sp|H2|H2_HUMAN\nMSAMPLEK\n>sp|Y1|Y1_YEAST\nMPEPTIDEKYEASTK\n', encoding='utf-8')
        self.registry = SpeciesRegistry(str(root / 'human.fasta'), str(root / 'both.fasta'))

    def tearDown(self):
        self.temp.cleanup()

    def test_one_way_rule(self):
        self.assertEqual(classify('H1', self.registry.target, self.registry.entrap), 'target')
        self.assertEqual(classify('Y1', self.registry.target, self.registry.entrap), 'entrapment')
        self.assertEqual(classify('H1;Y1', self.registry.target, self.registry.entrap), 'shared')
        self.assertEqual(classify('X9', self.registry.target, self.registry.entrap), 'unknown')
        # A yeast-only peptide found in a human protein when I is read as L becomes I/L-shared; human ones are not checked.
        self.assertEqual(self.registry.labels(['Y1', 'Y1', 'H1'], ['PEPTIDEK', 'YEASTK', 'PEPTLDEK']), ['shared_IL', 'entrapment', 'target'])

    def test_run_names_and_acceptance(self):
        self.assertEqual(run_name('runs/20230929_yzc_Hela_1.raw'), '20230929_yzc_Hela_1')
        frame = pd.DataFrame({'Run': ['20230929_yzc_Hela_1.raw'] * 3, 'Modified.Sequence': ['A', 'B', 'C'], 'Precursor.Charge': [2, 3, 2],
                              'Q.Value': [0.01, 0.011, 0.001], 'Decoy': [0, 0, 1]})
        accepted = accept_rows(frame, targets_only=True)
        self.assertEqual(accepted['key'].tolist(), [('A', 2)])
        self.assertEqual(accepted['Run'].tolist(), ['Hela_1'])

    def test_counts_and_ratios(self):
        record = species_record(['target', 'target', 'entrapment', 'shared', 'shared_IL', 'unknown'])
        self.assertEqual((record['human'], record['yeast'], record['species_exclusive_denominator']), (2, 1, 3))
        self.assertAlmostEqual(record['observed_yeast_fraction'], 1 / 3)
        self.assertEqual(fraction_ratio(0.0, 0.0) != fraction_ratio(0.0, 0.0), True)
        self.assertEqual(fraction_ratio(0.1, 0.0), np.inf)
        self.assertEqual(collapse_labels(['target', 'entrapment']), 'shared')
        self.assertAlmostEqual(combined_factor(0.25), 5.0)


class SupportAndQuantity(unittest.TestCase):
    def test_support_bins(self):
        self.assertEqual([support_bin(k, 63) for k in [0, 12, 13, 31, 32, 63]], ['zero', 'low_positive', 'middle', 'middle', 'high', 'high'])
        self.assertEqual([original_bin(k, 10) for k in [0, 2, 4, 6, 7]], ['zero', 'positive_0_20', 'positive_20_40', 'positive_40_60', 'positive_60_100'])

    def test_terciles_use_each_cells_first_pass_cutoffs(self):
        first = pd.DataFrame({'Run': ['a'] * 3 + ['b'] * 3, 'Precursor.Quantity': [1., 2., 3., 10., 20., 30.]})
        selected = pd.DataFrame({'Run': ['a', 'a', 'b', 'b'], 'Precursor.Quantity': [1.5, 3.5, 5., np.nan]})
        terciles, cuts = quantity_terciles(first, selected)
        self.assertEqual(list(terciles), ['low', 'high', 'low', 'quantity_unavailable'])
        self.assertEqual(len(cuts), 2)

    def test_wilson(self):
        low, high = wilson(0, 100)
        self.assertAlmostEqual(low, 0.0, places=12)
        self.assertTrue(0.03 < high < 0.04)


class Intervals(unittest.TestCase):
    def test_percentile_interval_and_undefined_draws(self):
        result = bb.interval(np.arange(1, 2001, dtype=float))
        self.assertAlmostEqual(result['low95'], 50.975)
        self.assertAlmostEqual(result['high95'], 1950.025)
        self.assertEqual(bb.interval([1.0, np.nan])['interval_status'], 'not estimable: undefined draw')

    def test_ratio_handles_zero_reference(self):
        values = bb.ratio(np.array([1., 0.]), np.array([10., 10.]), np.array([0., 1.]), np.array([10., 10.]))
        self.assertTrue(np.isposinf(values[0]))
        self.assertEqual(values[1], 0.0)

    def test_multiplicities_sum_to_units(self):
        weights = bb.multiplicities(np.random.default_rng(1), 7, 5)
        self.assertTrue((weights.sum(axis=1) == 7).all())

    def test_derived_classes_add_up(self):
        base = {'retained': (1, 10), 'second_pass_only': (2, 20), 'lost_library_present': (3, 30), 'lost_library_excluded': (4, 40)}
        out = derive_counts(base)
        self.assertEqual(out['lost'], (7, 70))
        self.assertEqual(out['first_pass'], (8, 80))
        self.assertEqual(out['final'], (3, 30))

    def test_count_model_rate_ratio(self):
        point, low, high = rate_ratio_ci(10, 1000, 10, 1000)
        self.assertEqual(point, 1.0)
        self.assertTrue(low < 1 < high)


class Mapping(unittest.TestCase):
    def test_protein_groups_are_joined_when_they_share_a_precursor(self):
        keys = [('A', 2), ('B', 2), ('C', 2)]
        ann = {('A', 2): {'stripped': 'A', 'groups': ['G1', 'G2']}, ('B', 2): {'stripped': 'B', 'groups': ['G2']},
               ('C', 2): {'stripped': 'C', 'groups': []}}
        units = group_units(keys, ann, 'protein_group')
        self.assertEqual(units[('A', 2)], units[('B', 2)])
        self.assertTrue(units[('C', 2)].startswith('unannotated:'))

    def test_natural_order(self):
        self.assertEqual(sorted(['Hela_10', 'Hela_2', 'Hela_1'], key=natural), ['Hela_1', 'Hela_2', 'Hela_10'])


if __name__ == '__main__':
    unittest.main()
