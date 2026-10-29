"""Figure 2 data: every plotted value comes from the class counts and intervals, with whole passes summed from counts."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from _paths import REPO  # noqa: F401
from mbraudit.figures import figure2

# (group, class, A, P) of the 54-run analysis and of the human-only runs analysed on their own (Figure 2 of the article)
MIXED = [('H', 'R', 79, 17523), ('H', 'T', 220, 6488), ('H', 'L', 22, 3465), ('Y', 'R', 1809, 13598), ('Y', 'T', 1430, 3490),
         ('Y', 'L', 181, 2334), ('E', 'R', 1348, 16717), ('E', 'T', 1196, 3662), ('E', 'L', 106, 2340)]
ALONE = [('H', 'R', 79, 17902), ('H', 'T', 39, 7661), ('H', 'L', 29, 3097)]


def counts(rows, pooled):
    out = []
    for group, cls, a, p in rows:
        out.append({'endpoint': 'T1', 'mask': 'F1', 'scope': 'organism', 'group': group, 'class': cls, 'tercile': 'all', 'A': a, 'P': p})
    if pooled:
        for cls in 'RTL':
            a = sum(r[2] for r in rows if r[1] == cls)
            p = sum(r[3] for r in rows if r[1] == cls)
            out.append({'endpoint': 'T1', 'mask': 'F1', 'scope': 'pooled', 'group': 'all', 'class': cls, 'tercile': 'all', 'A': a, 'P': p})
    # rows of another filter and of run scope must be ignored
    out.append({'endpoint': 'T1', 'mask': 'F2', 'scope': 'organism', 'group': 'H', 'class': 'T', 'tercile': 'all', 'A': 1, 'P': 1})
    out.append({'endpoint': 'T1', 'mask': 'F1', 'scope': 'run', 'group': 'H', 'class': 'T', 'tercile': 'all', 'A': 2, 'P': 3})
    frame = pd.DataFrame(out)
    frame['n'] = frame.A + frame.P
    frame['absent_fraction'] = frame.A / frame.n
    frame.insert(0, 'endpoint_id', range(len(frame)))
    return frame


class PanelData(unittest.TestCase):
    def test_plotted_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            mixed, alone = Path(temporary) / 'mixed', Path(temporary) / 'alone'
            mixed.mkdir()
            alone.mkdir()
            table = counts(MIXED, pooled=True)
            table.to_csv(mixed / 'endpoint_counts.tsv', sep='\t', index=False)
            intervals = pd.DataFrame({'endpoint_id': list(table.endpoint_id) * 2, 'scheme': ['two_way'] * len(table) + ['recipient'] * len(table),
                                      'low95': [i / 1000 for i in table.endpoint_id] + [9.0] * len(table),
                                      'high95': [i / 100 for i in table.endpoint_id] + [9.0] * len(table)})
            intervals.to_csv(mixed / 'fraction_intervals.tsv', sep='\t', index=False)
            counts(ALONE, pooled=False).to_csv(alone / 'endpoint_counts.tsv', sep='\t', index=False)
            frame = figure2.panel_data(mixed, alone)
        a = frame[frame.panel == 'A'].set_index(['x', 'class'])
        self.assertEqual(len(a), 9)  # pooled rows are not plotted
        added = a.loc[('Human-only runs', 'Added')]
        self.assertEqual((added.A, added.P, added.n), (220, 6488, 6708))
        self.assertEqual(added.percent, 100 * (220 / 6708))
        self.assertEqual((added.low95_percent, added.high95_percent), (100 * 0.001, 100 * 0.01))  # two-way limits of endpoint 1
        self.assertAlmostEqual(a.loc[('E. coli-only runs', 'Lost')].percent, 100 * 106 / 2446, places=12)
        b = frame[frame.panel == 'B']
        self.assertEqual(list(dict.fromkeys(b['class'])), ['First pass', 'Retained', 'Added', 'Lost', 'Final result'])
        b = b.set_index(['class', 'analysis'])
        first = b.loc[('First pass', 'with the other 48 runs')]
        self.assertEqual((first.A, first.n), (101, 21089))  # retained plus lost counts, not averaged fractions
        self.assertEqual(round(first.percent, 3), 0.479)
        final = b.loc[('Final result', 'on their own')]
        self.assertEqual((final.A, final.P, final.n), (118, 25563, 25681))
        self.assertEqual(round(final.percent, 3), 0.459)
        self.assertEqual(b.loc[('Added', 'on their own')].A, 39)
        self.assertTrue(b.low95_percent.isna().all())


if __name__ == '__main__':
    unittest.main()
