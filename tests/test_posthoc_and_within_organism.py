"""Hand-calculated cases for the post hoc fold-change analyses and the within-organism accounting
(adapted from the earlier release's tests; the expected values are unchanged)."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from _paths import REPO  # noqa: F401
from mbraudit import pipeline
from mbraudit.pg250 import posthoc_fc as c
from mbraudit.pg250 import within_organism as s


def posthoc_fixture(root):
    root.mkdir(parents=True)
    samples = pd.DataFrame([{'sample': 'SYN_' + g + '_' + str(i), 'kind': 'mixture', 'group': g} for g in ['5H4Y1E', '5H1Y4E'] for i in range(6)])
    names = [('SYN_' + sp + '_' + str(i), 2) for sp in ['E', 'Y'] for i in range(6)]
    keys = pd.DataFrame([(i, *key) for i, key in enumerate(names)], columns=['key', 'Modified.Sequence', 'Precursor.Charge'])
    rows, events = [], []
    for key, seq, charge in keys.itertuples(index=False, name=None):
        sp = seq.split('_')[1]
        i = int(seq.split('_')[2])
        sign = 1 if sp == 'Y' else -1
        scale = 2. ** (-12 + 4 * i)
        definitions = {'D': (2, 2, 2.), 'R': (2, 2, 2.), 'T': (2, 2, 1.), 'F': (4, 4, 1.5)}
        if i == 1:
            definitions.update(T=(0, 0, np.nan), F=(2, 2, 2.))
        if i == 2:
            definitions.update(D=(1, 1, np.nan), R=(1, 1, np.nan), F=(3, 3, 1.))
        if i == 3:
            definitions.update(T=(1, 0, np.nan), F=(3, 2, 2.))
        for mask in c.MASKS:
            for quantity in c.QUANTITIES:
                for context, (nh, nl, fc) in definitions.items():
                    observed = sign * fc
                    rows.append(dict(mask=mask, species=sp, context=context, quantity=quantity, key=key, batch='SYN_DATE', n_high=nh, n_low=nl,
                                     eligible=nh >= 2 and nl >= 2, observed_log2fc=observed, expected_log2fc=sign * 2., error=observed - sign * 2.,
                                     coverage_population='own_context'))
            for g in ['5H4Y1E', '5H1Y4E']:
                high = g == '5H4Y1E'
                r_count = 1 if i == 2 else 2
                t_count = 0 if i == 1 else (int(high) if i == 3 else 2)
                for position in range(r_count + t_count):
                    ratio = 4 if position < r_count else 2
                    quantity = scale * (ratio if (high and sp == 'Y') or (not high and sp == 'E') else 1)
                    events.append({'sample': 'SYN_' + g + '_' + str(position), 'Modified.Sequence': seq, 'Precursor.Charge': charge, 'species': sp,
                                   'mask': mask, 'stage': 'final', 'Precursor.Normalised': quantity})
    frames = {'fc': pd.DataFrame(rows), 'keys': keys, 'events': pd.DataFrame(events), 'samples': samples}
    studies = {}
    for study, p3 in [('SYN_A', True), ('SYN_B', False)]:
        inputs = {}
        for role in (frames if p3 else ['fc', 'keys']):
            path = root / ('SYN_' + study + '_' + role + '.tsv')
            frames[role].to_csv(path, sep='\t', index=False, lineterminator='\n', na_rep='NA')
            inputs[role] = path
        studies[study] = {'P3': p3, 'batch_date': 'SYN_DATE', 'inputs': inputs}
    return studies


class Posthoc(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        root = Path(cls.temp.name)
        pipeline.posthoc(posthoc_fixture(root / 'inputs'), root / 'results', synthetic=True)
        cls.summary = pd.read_csv(root / 'results' / 'summary.tsv', sep='\t')

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def row(self, endpoint, metric, species='Y', **fields):
        f = self.summary
        mask = (f.batch == 'SYN_A') & (f['mask'] == 'F1') & (f.quantity == 'Precursor.Normalised') & (f.species == species) & (f.endpoint == endpoint) & (f.metric == metric)
        for key, value in fields.items():
            mask &= f[key].eq(value)
        return f[mask].iloc[0]

    def test_paired_population_and_error_change(self):
        self.assertEqual(self.row('P1', 'median_delta').eligible_keys, 4)
        self.assertEqual(self.row('P1', 'median_delta').point, -.5)
        self.assertEqual(self.row('P1', 'MAE_change').point, .375)
        self.assertEqual(self.row('P1', 'fraction_toward_zero').point, .75)

    def test_abundance_terciles_are_fixed_on_the_union_before_the_class_split(self):
        for tier, r, t, shared in [('low', 2, 1, 1), ('middle', 1, 1, 0), ('high', 2, 2, 2)]:
            row = self.row('P3', 'median_bias_T_minus_R', tercile=tier)
            self.assertEqual((row.union_keys, row.R_keys, row.T_keys, row.shared_keys), (2, r, t, shared))
            self.assertEqual(row.point, -1.)


class WithinOrganism(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / 'inputs'
        self.root.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def fixture(self, organism='H'):
        names = ['SYN_' + organism + '_' + str(i) for i in range(6)]
        design = pd.DataFrame({'sample': names, 'kind': 'single', 'group': organism, 'filename_date': 'SYN_DATE', 'present_species': organism})
        absent = next(x for x in s.BITS_SPECIES if x != organism)
        annotation = pd.DataFrame([{s.KEY[0]: key, s.KEY[1]: 2, 'species': sp, 'Stripped.Sequence': key, 'Protein.Ids': 'SYN_PROTEIN_' + key}
                                   for key, sp in [('SYN_A', organism), ('SYN_B', absent), ('SYN_T', absent), ('SYN_U', 'unknown')]])

        def row(sample, key, decoy=0, quantity=10):
            return {'sample': sample, s.KEY[0]: key, s.KEY[1]: 2, 'Decoy': decoy, 'Precursor.Quantity': quantity, 'Precursor.Normalised': quantity,
                    **{field: .001 for field in s.MASKS['F2']}}
        first = [row(n, 'SYN_A') for n in names[:5]] + [row(names[i], 'SYN_B') for i in [0, 1]] + [row(names[0], 'SYN_A', 1, 1000)]
        final = [row(names[i], 'SYN_A') for i in [0, 2, 3, 4]] + [row(names[0], 'SYN_B'), row(names[0], 'SYN_T', quantity=20), row(names[2], 'SYN_U', quantity=0)]
        empirical = pd.DataFrame([{s.KEY[0]: 'SYN_A', s.KEY[1]: 2, 'Decoy': 0}])
        exports = pd.DataFrame([{'sample': name, 'stage': stage, 'status': 'complete'} for name in names for stage in ['first', 'final']])
        files = {'first': self.root / 'report-first-pass.parquet', 'final': self.root / 'report.parquet', 'empirical': self.root / 'report-lib.parquet',
                 'predicted_species': self.root / 'SYN_species.parquet', 'design': self.root / 'SYN_design.tsv', 'exports': self.root / 'SYN_exports.tsv'}
        for role, frame in [('first', pd.DataFrame(first)), ('final', pd.DataFrame(final)), ('empirical', empirical), ('predicted_species', annotation)]:
            frame.to_parquet(files[role], index=False)
        design.to_csv(files['design'], sep='\t', index=False)
        exports.to_csv(files['exports'], sep='\t', index=False)
        return files, design

    def test_six_runs_decoys_and_hand_counted_classes(self):
        files, design = self.fixture()
        out = Path(self.temp.name) / 'results'
        pipeline.within_organism(files, out, 'H', synthetic=True)
        summary = pd.read_csv(out / 'endpoint_counts.tsv', sep='\t')
        rows = summary[(summary['mask'] == 'F1') & (summary.endpoint == 'T1') & (summary.scope == 'organism')].set_index('class')
        self.assertEqual(rows.events.to_dict(), {'D': 7, 'F': 7, 'R': 5, 'T': 2, 'L': 2, 'L_in_library': 1, 'L_outside_library': 1})
        self.assertEqual((rows.loc['D', 'A'], rows.loc['D', 'P']), (2, 5))
        self.assertEqual((rows.loc['T', 'A'], rows.loc['T', 'P'], rows.loc['T', 'excluded_unknown']), (1, 0, 1))
        self.assertEqual(len(summary[summary.group.eq(design['sample'].iloc[5])]), 81)
        self.assertEqual(set(summary.scope), {'run', 'organism'})


if __name__ == '__main__':
    unittest.main()
