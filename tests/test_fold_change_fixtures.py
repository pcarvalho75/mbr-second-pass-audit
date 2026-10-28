"""Completeness selection, fold changes and whole-pass arithmetic on the released synthetic fixtures.

The two fixtures (accepted events, designed samples, run status) and their expected output hashes were released
with the earlier code version; the recomputed tables must match them byte for byte. The hashes were refreshed once,
when the filter labels were renamed to those of the article (F3 -> F1, F4 -> F2, F0 -> F3); the tables are otherwise
those of the earlier version (same rows and values; fc_summary.tsv and t3_summary.tsv are sorted by the label).
"""
import csv
import json
import tempfile
import unittest
from pathlib import Path

from _paths import FIXTURES
from mbraudit import pipeline
from mbraudit.pg250.guard import sha, verify_inputs
from mbraudit.pg250.summary import accepted_passes, added_count_comparison
from mbraudit.singlecell.scales import final_single_cell_bound, scale_point


class FixtureOutputs(unittest.TestCase):
    def run_fixture(self, study):
        folder = FIXTURES / 'accepted-events' / study
        inputs = {n: folder / (n + '.tsv') for n in ['events', 'samples', 'status']}
        names = json.loads((folder / 'sample-names.json').read_text())['sample_names']
        expected = json.loads((folder / 'expected-outputs.json').read_text())['outputs']
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / 'results'
            pipeline.fold_changes(inputs, out, study, mode='synthetic', sample_names=names)
            for name, item in expected.items():
                self.assertEqual(sha(out / name), item['sha256'], name)

    def test_primary_tables_match_all_eight_reference_hashes(self):
        self.run_fixture('primary')

    def test_replication_tables_match_all_four_reference_hashes(self):
        self.run_fixture('replication')

    def test_stale_input_is_refused(self):
        path = FIXTURES / 'accepted-events' / 'primary' / 'events.tsv'
        with self.assertRaisesRegex(ValueError, 'Input identity differs'):
            verify_inputs({'events': path}, {'events': '0' * 64})


class Arithmetic(unittest.TestCase):
    def test_all_fifteen_upper_bounds_keep_their_precision(self):
        with (FIXTURES / 'search-space-scales.tsv').open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 15)
        for row in rows:
            y, h = map(int, row['r_exact'].split('/'))
            fraction, percent, display = scale_point(row['full_precision_f'], h, y, row['configuration'], row['class'])
            self.assertEqual((fraction, percent, display),
                             (row['search_space_scale_fraction'], row['search_space_scale_percent'], row['Search-space scale']))

    def test_whole_final_result_counts_retained_and_added_once(self):
        rows = [{'filter': 'F1', 'class': 'retained', 'yeast': '164', 'species_exclusive_denominator': '246201'},
                {'filter': 'F1', 'class': 'second_pass_only', 'yeast': '218', 'species_exclusive_denominator': '100239'}]
        result = final_single_cell_bound(rows, 1365729, 374286)
        self.assertEqual((result['yeast'], result['denominator']), (382, 346440))
        self.assertEqual((result['lower_percent'], result['estimated_upper_percent']), ('0.110', '0.513'))

    def test_table4_sums_counts_rather_than_averaging_fractions(self):
        def rows(values):
            return [{'endpoint': 'T1', 'mask': 'F1', 'scope': 'organism', 'group': 'H', 'tercile': 'all',
                     'class': cls, 'A': str(a), 'P': str(p)} for cls, a, p in values]
        mixed = rows([('R', 79, 17523), ('T', 220, 6488), ('L', 22, 3465)])
        own = rows([('R', 79, 17902), ('T', 39, 7661), ('L', 29, 3097)])
        self.assertEqual([(r['A'], r['P'], r['n'], r['percent']) for r in accepted_passes(mixed, 'H')],
                         [(101, 20988, 21089, '0.479'), (299, 24011, 24310, '1.230')])
        self.assertEqual([(r['A'], r['P'], r['n'], r['percent']) for r in accepted_passes(own, 'H')],
                         [(108, 20999, 21107, '0.512'), (118, 25563, 25681, '0.459')])
        change = added_count_comparison(mixed, own)
        self.assertEqual((change['wrong_added_mixed'], change['wrong_added_own']), (220, 39))
        self.assertEqual(round(float(change['wrong_added_count_reduction_percent'])), 82)


if __name__ == '__main__':
    unittest.main()
