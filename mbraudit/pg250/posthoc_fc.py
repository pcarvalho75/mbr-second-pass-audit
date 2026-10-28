"""Post hoc fold-change analyses P1 (paired), P2 (call consequences) and P3 (abundance terciles) (carried from validation_v5/competitive_fc.py)."""
from pathlib import Path
import argparse, hashlib, importlib.metadata, json, warnings
import numpy as np
import pandas as pd
from mbraudit.pg250.fold_change_replication import interval, SEED
DRAWS = 2000
QUANTITIES = ('Precursor.Normalised', 'Precursor.Quantity')
MASKS = ('F3', 'F0', 'F4')
_SEAL = object()

def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')

def bootstrap(arrays, statistic):
    n = len(arrays[0])
    rng = np.random.default_rng(SEED)
    values = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        point = float(statistic(*[a[None, :] for a in arrays])[0]) if n else np.nan
        for start in range(0, DRAWS, 16):
            indices = rng.integers(0, n, size=(min(16, DRAWS - start), n)) if n else np.empty((min(16, DRAWS - start), 0), int)
            values.extend(statistic(*[a[indices] for a in arrays]))
    return (point, interval(values))

def numeric(series):
    return series.map(lambda x: np.nan if x == 'NA' else float(x))

def analyse(guard):
    if not isinstance(guard, dict) or guard.get('seal') is not _SEAL:
        raise ValueError('Verified analysis guard required')
    summaries = []
    members = []
    cells = []
    ledgers = []
    for study, scope in guard['config']['studies'].items():
        frames = {}
        for role in scope['inputs']:
            alias = study + '.' + role
            raw = pd.read_csv(guard['inputs'][alias], sep='\t', dtype=str, keep_default_na=False)
            if any(('\n' in str(x) or '\r' in str(x) for x in raw.to_numpy().flat)):
                raise ValueError('Multiline TSV cell prevents physical-line ledger')
            sha = scope['inputs'][role]['sha256']
            for line, row in enumerate(raw.to_dict('records'), 2):
                cells.append({'source_ref': alias + ':' + str(line), 'source_sha256': sha, 'line': line, 'columns': json.dumps(row, sort_keys=True)})
            raw['_ref'] = [alias + ':' + str(i) for i in range(2, len(raw) + 2)]
            frames[role] = raw
        if guard['synthetic'] and (not frames['keys']['Modified.Sequence'].str.startswith('SYN_').all() or any((not f['sample'].str.startswith('SYN_').all() for f in frames.values() if 'sample' in f))):
            raise ValueError('Non-SYN identity')
        fc = frames['fc'].copy()
        keys = frames['keys'].copy()
        for c in ['key', 'n_high', 'n_low', 'observed_log2fc', 'expected_log2fc', 'error']:
            fc[c] = numeric(fc[c])
        keys['key'] = numeric(keys.key)
        keys['Precursor.Charge'] = numeric(keys['Precursor.Charge'])
        if keys.key.duplicated().any() or keys.duplicated(['Modified.Sequence', 'Precursor.Charge']).any() or (keys.key % 1).ne(0).any():
            raise ValueError('Duplicate/invalid canonical key')
        if not set(fc.key) <= set(keys.key) or fc.duplicated(['mask', 'species', 'quantity', 'context', 'key']).any():
            raise ValueError('Duplicate/unknown FC key')
        if not set(fc['mask']) <= set(MASKS) or not set(fc.species) <= {'Y', 'E'} or (not set(fc.context) <= set('DFRT')) or (not set(fc.quantity) <= set(QUANTITIES)):
            raise ValueError('FC category differs')
        if set(fc.eligible) - {'True', 'False'} or (fc[['key', 'n_high', 'n_low']] % 1).ne(0).any().any() or fc[['n_high', 'n_low']].lt(0).any().any() or fc[['n_high', 'n_low']].gt(6).any().any():
            raise ValueError('Invalid eligibility/count')
        if set(fc.batch) != {scope['batch_date']}:
            raise ValueError('Batch label differs')
        fc['eligible'] = fc.eligible.map({'True': True, 'False': False})
        expected = fc.species.map({'Y': 2.0, 'E': -2.0})
        if not fc.expected_log2fc.eq(expected).all() or not fc.eligible.eq((fc.n_high >= 2) & (fc.n_low >= 2)).all():
            raise ValueError('Eligibility/truth differs')
        if not np.isfinite(fc.loc[fc.eligible, 'observed_log2fc']).all() or fc.loc[~fc.eligible, 'observed_log2fc'].notna().any():
            raise ValueError('Invalid eligible quantity')
        if not fc.loc[fc.eligible, 'error'].eq(fc.loc[fc.eligible, 'observed_log2fc'] - expected[fc.eligible]).all():
            raise ValueError('Frozen error inconsistency')
        keyrefs = keys.set_index('key')['_ref'].to_dict()
        keynames = keys.set_index('key')['Modified.Sequence'].to_dict()
        if scope['P3']:
            samples = frames['samples']
            events = frames['events'].copy()
            axis = samples[(samples.kind == 'mixture') & samples['group'].isin(['5H4Y1E', '5H1Y4E'])]
            if samples['sample'].duplicated().any() or axis.groupby('group').size().to_dict() != {'5H1Y4E': 6, '5H4Y1E': 6}:
                raise ValueError('Extreme design differs')
            events['Precursor.Charge'] = numeric(events['Precursor.Charge'])
            events['Precursor.Normalised'] = numeric(events['Precursor.Normalised'])
            if events.duplicated(['sample', 'mask', 'stage', 'Modified.Sequence', 'Precursor.Charge']).any():
                raise ValueError('Duplicate accepted event')
            events = events.merge(keys[['key', 'Modified.Sequence', 'Precursor.Charge']], on=['Modified.Sequence', 'Precursor.Charge'], how='left', validate='many_to_one')
            if events.key.isna().any() or not set(events['sample']) <= set(samples['sample']):
                raise ValueError('Unmapped accepted event')
            events = events[(events.stage == 'final') & events['sample'].isin(axis['sample']) & events['Precursor.Normalised'].gt(0) & np.isfinite(events['Precursor.Normalised'])].copy()
            events['log_abundance'] = np.log2(events['Precursor.Normalised'])
        for mask in MASKS:
            for sp in ['Y', 'E']:
                for quantity in QUANTITIES:
                    g = fc[(fc['mask'] == mask) & (fc.species == sp) & (fc.quantity == quantity)]
                    universe = sorted(g.key.unique())
                    truth = 2.0 if sp == 'Y' else -2.0
                    piv = {c: g[g.context == c].set_index('key').reindex(universe) for c in 'DFRT'}
                    if len(g) and any((p['context'].isna().any() for p in piv.values())):
                        raise ValueError('Incomplete frozen context universe')
                    values = {c: p.observed_log2fc.to_numpy(float) for c, p in piv.items()}
                    elig = {c: p.eligible.fillna(False).to_numpy(bool) for c, p in piv.items()}
                    base = {'batch': study, 'mask': mask, 'species': sp, 'quantity': quantity, 'post_hoc': True, 'display': 'primary' if mask == 'F3' and quantity == QUANTITIES[0] else 'sensitivity'}

                    def population(kind, selection, extra=None):
                        name = '|'.join([study, mask, sp, quantity, kind])
                        for i, k in enumerate(universe):
                            row = {**base, 'population': name, 'key': int(k), 'Modified.Sequence': keynames[k], 'selected': bool(selection[i]), 'key_ref': keyrefs[k], **{c + '_ref': piv[c].iloc[i]['_ref'] for c in 'DFRT'}}
                            if extra:
                                row.update({field: array[i] for field, array in extra.items()})
                            members.append(row)
                        return name

                    def emit(endpoint, metric, name, arrays, statistic, counts):
                        point, ci = bootstrap(arrays, statistic)
                        row = {**base, 'endpoint': endpoint, 'metric': metric, 'population': name, 'universe_keys': len(universe), **counts, 'point': point, 'ci_low': ci[0], 'ci_high': ci[1], 'undefined_draws': int(ci[2]), 'infinite_draws': int(ci[3]), 'status': 'estimable' if np.isfinite(point) else 'not_estimable'}
                        summaries.append(row)
                        ledgers.append({'summary_row': len(summaries) + 1, 'population': name, 'metric': metric, 'formula': metric, 'source_table': 'membership.tsv', 'source_selection': 'all rows with this population; selected and component flags fix masks', 'universe_order': 'key ascending', 'bootstrap': 'PCG64 reset 20261018; 2000 draws; batches16; full universe', 'interval': 'linear .025/.975; unavailable if any undefined/infinite draw'})
                    mean = lambda a: np.nanmean(a, axis=1)
                    median = lambda a: np.nanmedian(a, axis=1)
                    selected = elig['R'] & elig['F'] & ((piv['T'].n_high.fillna(0) + piv['T'].n_low.fillna(0)).to_numpy() > 0)
                    delta = values['F'] - values['R']
                    change = np.abs(values['F'] - truth) - np.abs(values['R'] - truth)
                    toward = np.abs(values['F']) < np.abs(values['R'])
                    ties = np.abs(values['F']) == np.abs(values['R'])
                    name = population('P1', selected, {'FC_R': values['R'], 'FC_F': values['F'], 'delta': delta, 'absolute_error_change': change, 'toward_zero': toward, 'tie': ties, 'away_from_zero': ~toward & ~ties & np.isfinite(values['R']) & np.isfinite(values['F'])})
                    for metric, array, fn in [('median_delta', delta, median), ('fraction_toward_zero', toward, mean), ('MAE_change', change, mean), ('MAE_R', np.abs(values['R'] - truth), mean), ('MAE_F', np.abs(values['F'] - truth), mean), ('fraction_tie', ties, mean), ('fraction_away_from_zero', ~toward & ~ties, mean)]:
                        emit('P1', metric, name, [np.where(selected, array, np.nan)], fn, {'eligible_keys': int(selected.sum())})
                    for c in 'DFRT':
                        name = population('P2_' + c, elig[c], {'context': np.full(len(universe), c), 'FC': values[c], 'missed_twofold': np.abs(values[c]) < 1, 'wrong_sign': values[c] * truth < 0, 'within_half': np.abs(values[c] - truth) <= 0.5})
                        for metric, array in [('missed_twofold', np.abs(values[c]) < 1), ('wrong_sign', values[c] * truth < 0), ('within_half', np.abs(values[c] - truth) <= 0.5)]:
                            emit('P2', metric, name, [np.where(elig[c], array, np.nan)], mean, {'context': c, 'eligible_keys': int(elig[c].sum()), 'numerator': int((array & elig[c]).sum())})
                    if not scope['P3'] or quantity != QUANTITIES[0]:
                        continue
                    union = elig['R'] | elig['T']
                    part = events[(events['mask'] == mask) & (events.species == sp) & events.key.isin(np.array(universe)[union])]
                    abundance = part.groupby('key').log_abundance.median().reindex(universe)
                    if abundance[union].isna().any():
                        raise ValueError('Missing union abundance')
                    ranks = (abundance[union].rank(method='average') - 0.5) / int(union.sum()) if union.any() else pd.Series(dtype=float)
                    bins = pd.Series('not_in_union', index=universe, dtype='str')
                    bins.loc[ranks.index] = np.where(ranks <= 1 / 3, 'low', np.where(ranks <= 2 / 3, 'middle', 'high'))
                    refs = part.groupby('key')['_ref'].agg(lambda x: ','.join(x)).reindex(universe).fillna('').to_numpy()
                    for tier in ['low', 'middle', 'high']:
                        keep = bins.eq(tier).to_numpy()
                        left = keep & elig['R']
                        right = keep & elig['T']
                        name = population('P3_' + tier, keep, {'abundance': abundance.to_numpy(), 'tercile': bins.to_numpy(), 'eligible_R': elig['R'], 'eligible_T': elig['T'], 'abundance_refs': refs})
                        a = np.where(left, values['R'] - truth, np.nan)
                        b = np.where(right, values['T'] - truth, np.nan)
                        count = {'tercile': tier, 'union_keys': int(keep.sum()), 'R_keys': int(left.sum()), 'T_keys': int(right.sum()), 'shared_keys': int((left & right).sum())}
                        for metric, arrays, fn in [('median_bias_R', [a], median), ('median_bias_T', [b], median), ('median_bias_T_minus_R', [a, b], lambda x, y: np.nanmedian(y, axis=1) - np.nanmedian(x, axis=1))]:
                            emit('P3', metric, name, arrays, fn, count)
    return {'summary': pd.DataFrame(summaries), 'membership': pd.DataFrame(members), 'source_cells': pd.DataFrame(cells), 'summary_ledger': pd.DataFrame(ledgers)}
