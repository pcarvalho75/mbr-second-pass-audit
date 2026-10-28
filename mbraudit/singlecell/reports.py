"""Read the first-pass and final DIA-NN reports of one search and label species.

Carried from analysis/revision_reports.py, with natural (analysis/cluster_v62.py), with_prec
(analysis/restoration_math.py), annotations and group_units (analysis/revision_restoration.py).
report_dir and Reports.pin are plumbing.
"""
from pathlib import Path
import hashlib
import re
import numpy as np
import pandas as pd
from mbraudit.species import accept_rows, SpeciesRegistry, ENTRAP_RATIO, poisson_combined_bounds, run_name, PXD_RUN_PREFIX
from mbraudit.filters import filter_mask
from mbraudit.singlecell.design import H6_RUNS, M6, H16, M16
from mbraudit import config


def report_dir(alias):
    """Folder of one search: all_cells_64 is the active configuration, other aliases are pool searches."""
    return config.single_cell_report_dir(alias)


KEY = ['Modified.Sequence', 'Precursor.Charge']


QCOLS = ['Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value', 'Lib.PG.Q.Value']


COLS = ['Run', *KEY, 'Decoy', 'Stripped.Sequence', 'Protein.Ids', 'Protein.Group', *QCOLS]


SEED = 20261018


DRAWS = 2000


class Reports:

    def __init__(self):
        self.registry = SpeciesRegistry()
        self.label_cache = {}
        self.cache = {}
        self.pins = {}

    def pin(self, path, alias):
        """Record the SHA-256 of a scientific input (plumbing; replaces the old archive writer)."""
        path = Path(path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        self.pins[alias] = {'sha256': digest}
        return path

    def label(self, data):
        pairs = list(zip(data['Protein.Ids'].fillna(''), data['Stripped.Sequence'].fillna('')))
        missing = sorted(set(pairs) - set(self.label_cache))
        if missing:
            self.label_cache.update(zip(missing, self.registry.labels([p[0] for p in missing], [p[1] for p in missing])))
        data['species'] = [self.label_cache[p] for p in pairs]
        return data

    def load(self, alias):
        if alias in self.cache:
            return self.cache[alias]
        folder = report_dir(alias)
        raw = {}
        for stage, file in [('first', 'report-first-pass.parquet'), ('final', 'report.parquet')]:
            path = self.pin(folder / file, alias + ':' + stage)
            frame = pd.read_parquet(path, columns=COLS)
            targets = frame[frame.Decoy.eq(0)].copy()
            targets['Run'] = targets.Run.map(run_name).str.replace(PXD_RUN_PREFIX, '', regex=False)
            targets['key'] = list(zip(targets[KEY[0]], targets[KEY[1]].astype(int)))
            raw[stage + '_raw'] = targets
            accepted = accept_rows(frame, targets_only=True)
            if accepted.duplicated(['Run', *KEY]).any():
                raise ValueError('Duplicate accepted precursor-cell event')
            raw[stage] = self.label(accepted)
        first = raw['first']
        own = {cell: set(g.key) for cell, g in first.groupby('Run')}
        counts = {}
        for keys in own.values():
            for key in keys:
                counts[key] = counts.get(key, 0) + 1
        final = raw['final']
        final['direct'] = [key in own.get(cell, set()) for cell, key in zip(final.Run, final.key)]
        final['support'] = [counts.get(key, 0) - int(key in own.get(cell, set())) for cell, key in zip(final.Run, final.key)]
        final['pool'] = alias
        raw['final'] = final
        raw['own'] = own
        self.cache[alias] = raw
        return raw

    def arm(self, aliases, filter_name='F0', first=False):
        parts = []
        for alias in aliases:
            data = self.load(alias)
            frame = data['first' if first else 'final'].copy()
            frame['pool'] = alias
            if first:
                frame['direct'] = True
                frame['support'] = 0
            else:
                frame = frame.loc[filter_mask(frame[QCOLS].fillna(2.0), filter_name)].copy()
            parts.append(frame)
        result = pd.concat(parts, ignore_index=True)
        assert not result.duplicated(['Run', *KEY]).any()
        return result

    def library(self):
        path = self.pin(report_dir('all_cells_64') / 'report-lib.parquet', 'all_cells_64:empirical_library')
        return pd.read_parquet(path, columns=[*KEY, 'Decoy', 'Q.Value'])


def transfers(frame):
    return frame.loc[~frame.direct].copy()


def count_record(frame):
    e = int(frame.species.eq('entrapment').sum())
    t = int(frame.species.eq('target').sum())
    n = t + e
    point, lo, hi = poisson_combined_bounds(e, n) if n else (np.nan, np.nan, np.nan)
    return {'events': len(frame), 'human': t, 'yeast': e, 'error_denominator': n, 'observed_yeast_fraction': e / n if n else np.nan, 'estimated_combined_fdp': point, 'poisson_low95': lo, 'poisson_high95': hi}


def natural(run: str):
    return [int(x) if x.isdigit() else x for x in re.split('(\\d+)', run)]


def with_prec(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out['prec'] = list(zip(out['Modified.Sequence'].astype(str), out['Precursor.Charge'].astype(int)))
    out['rk'] = list(zip(out['Run'].astype(str), out['prec']))
    return out


def annotations(first, final, keys):
    both = pd.concat([first, final], ignore_index=True)
    result = {}
    for key, g in both[both.prec.isin(keys)].groupby('prec', sort=False):
        stripped = sorted(set(g['Stripped.Sequence'].dropna().astype(str)) - {''})
        if len(stripped) != 1:
            raise ValueError('Missing or conflicting stripped sequence')
        groups = sorted(set(g['Protein.Group'].dropna().astype(str)) - {''})
        result[key] = {'stripped': stripped[0], 'groups': groups}
    assert set(result) == set(keys)
    return result


def group_units(keys, ann, unit):
    if unit == 'stripped_sequence':
        return {key: ann[key]['stripped'] for key in keys}
    if unit != 'protein_group':
        raise ValueError('Unknown grouped split')
    parent = {g: g for row in ann.values() for g in row['groups']}

    def find(g):
        while parent[g] != g:
            parent[g] = parent[parent[g]]
            g = parent[g]
        return g
    for row in ann.values():
        groups = row['groups']
        for g in groups[1:]:
            a, b = (find(groups[0]), find(g))
            parent[max(a, b)] = min(a, b)
    result = {}
    for key in keys:
        groups = ann[key]['groups']
        result[key] = 'group:' + find(groups[0]) if groups else f'unannotated:{key[0]}|{key[1]}'
    return result
