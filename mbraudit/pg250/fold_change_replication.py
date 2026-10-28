"""Fold changes with a configurable mixture date, used for the 9 January replication (carried from validation_v5/validation_t3_fc.py)."""
from pathlib import Path
import argparse, hashlib, json, math, warnings
import numpy as np
import pandas as pd
from scipy import sparse
from mbraudit.pg250.guard import require_guard
SPECIES = ('H', 'Y', 'E', 'shared', 'shared_IL', 'unknown')
MASKS = ('F3', 'F1', 'F2')
STAGES = ('first', 'final')
KEY = ['Modified.Sequence', 'Precursor.Charge']
SEED = 20261018

def positive(series):
    x = pd.to_numeric(series, errors='coerce').to_numpy(float)
    return np.isfinite(x) & (x > 0)

def validate(events, samples, status, *, mixture_date='20250116'):
    cols = ['sample', *KEY, 'species', 'mask', 'stage', 'Precursor.Quantity', 'Precursor.Normalised']
    if not set(cols) <= set(events):
        raise ValueError('Missing accepted-event columns')
    e = events[cols].copy()
    s = samples.copy()
    if len(s) != 54 or s['sample'].duplicated().any():
        raise ValueError('Expected 54 unique designed runs')
    mix = s[s.kind == 'mixture']
    if len(mix) != 30 or set(mix.filename_date.astype(str)) != {mixture_date}:
        raise ValueError('C3 mixture batch')
    if sorted(mix.groupby('group').size()) != [6] * 5 or set(mix['group']) != {'5H1Y4E', '5H2Y3E', '4H2Y2E', '5H3Y2E', '5H4Y1E'}:
        raise ValueError('C3 ratio counts')
    singles = s[s.kind == 'single']
    if len(singles) != 18 or any((sum(singles.present_species == sp) != 6 for sp in ('H', 'Y', 'E'))) or sum(s.kind == 'blank') != 6:
        raise ValueError('C3 recipient/blank counts')
    if not set(e['sample']) <= set(s['sample']):
        raise ValueError('Undesigned sample')
    if not set(e.species) <= set(SPECIES) or not set(e['mask']) <= set(MASKS) or (not set(e.stage) <= set(STAGES)):
        raise ValueError('Invalid category')
    charges = pd.to_numeric(e['Precursor.Charge'], errors='raise')
    if not np.all(np.isfinite(charges) & (charges == np.floor(charges)) & (charges >= 2) & (charges <= 4)):
        raise ValueError('Invalid charge')
    e['Precursor.Charge'] = charges.astype(int)
    if e['Modified.Sequence'].isna().any() or (e['Modified.Sequence'].astype(str).str.len() == 0).any():
        raise ValueError('Missing key')
    if e.duplicated(['sample', 'mask', 'stage', *KEY]).any():
        raise ValueError('Duplicate accepted events')
    if (e.groupby(KEY).species.nunique() > 1).any():
        raise ValueError('Conflicting species mapping')
    for q in ['Precursor.Quantity', 'Precursor.Normalised']:
        e[q] = pd.to_numeric(e[q], errors='raise')
    if 'Decoy' in events and pd.to_numeric(events.Decoy, errors='raise').ne(0).any():
        raise ValueError('Decoy in target input')
    required = {(sample, st, mask) for sample in s['sample'] for st in STAGES for mask in MASKS}
    if set(status.columns) < {'sample', 'stage', 'mask', 'status'} or status.duplicated(['sample', 'stage', 'mask']).any():
        raise ValueError('Invalid run status')
    actual = set(status[['sample', 'stage', 'mask']].itertuples(index=False, name=None))
    if actual != required or not set(status.status) <= {'complete', 'unavailable'}:
        raise ValueError('Missing designed status axes')
    bad = status[status.status != 'complete']
    if len(e.merge(bad[['sample', 'stage', 'mask']], on=['sample', 'stage', 'mask'])):
        raise ValueError('Events on unavailable axis')
    keys = sorted(set(e[KEY].itertuples(index=False, name=None)))
    keymap = {k: i for i, k in enumerate(keys)}
    e['key'] = [keymap[k] for k in e[KEY].itertuples(index=False, name=None)]
    return (e, s, status, keys)

def transitions(e):
    """Join accepted sets before any completeness/quantity selection."""
    e = e.copy()
    e['transition'] = pd.Series(index=e.index, dtype='str')
    for mask in MASKS:
        d = set(e[(e['mask'] == mask) & (e.stage == 'first')][['sample', 'key']].itertuples(index=False, name=None))
        f = set(e[(e['mask'] == mask) & (e.stage == 'final')][['sample', 'key']].itertuples(index=False, name=None))
        ix = e['mask'] == mask
        e.loc[ix, 'transition'] = ['R' if (sample, key) in (f if stage == 'first' else d) else 'L' if stage == 'first' else 'T' for sample, key, stage in e.loc[ix, ['sample', 'key', 'stage']].itertuples(index=False, name=None)]
    return e

def interval(draws):
    a = np.asarray(draws, float)
    undefined = int(np.isnan(a).sum())
    infinite = int(np.isinf(a).sum())
    if undefined or infinite:
        return (np.nan, np.nan, undefined, infinite)
    return (*np.quantile(a, [0.025, 0.975], method='linear'), undefined, infinite)

def t3(e, s, status, keys, draws=2000, *, guard=None):
    require_guard(guard, producer=__file__, frames=(e, s, status))
    singles = s[s.kind == 'single'].sort_values('sample')
    runlist = list(singles['sample'])
    runidx = {r: i for i, r in enumerate(runlist)}
    subs = e[e['sample'].isin(runlist)].copy()
    keylist = sorted(subs.key.unique())
    keyidx = {k: i for i, k in enumerate(keylist)}
    records = []
    coverage = []
    membership = []
    contributions = []
    unavailable = []
    for mask in MASKS:
        for organism in ('H', 'Y', 'E'):
            axis = list(singles[singles.present_species == organism]['sample'])
            complete = {st: status[(status['mask'] == mask) & (status.stage == st) & status['sample'].isin(axis)].status.eq('complete').all() for st in STAGES}
            for scheme in ('stage_specific', 'common_first'):
                for stage in STAGES:
                    if not complete[stage] or not complete['first'] or (not complete['final']):
                        unavailable.append(dict(mask=mask, organism=organism, scheme=scheme, stage=stage, status='not_estimable', reason='required mask/stage unavailable for original partition or completeness'))
                        continue
                    base = subs[(subs['mask'] == mask) & (subs.stage == stage) & subs['sample'].isin(axis)]
                    selector = stage if scheme == 'stage_specific' else 'first'
                    obs = subs[(subs['mask'] == mask) & (subs.stage == selector) & subs['sample'].isin(axis)]
                    counts = obs[positive(obs['Precursor.Quantity'])].groupby('key')['sample'].nunique().to_dict()
                    for threshold in (0, 3, 5):
                        selected = set(base.key) if threshold == 0 else {k for k, n in counts.items() if n >= threshold}
                        if scheme == 'stage_specific' or stage == 'final':
                            for k in sorted(set(base.key) | set(counts)):
                                membership.append(dict(mask=mask, organism=organism, scheme=scheme, stage=stage, selector_stage=selector, threshold=threshold, key=k, observed_runs=counts.get(k, 0), designed_runs=6, selected=k in selected))
                        kept = base[base.key.isin(selected)]
                        classes = ('D', 'R', 'L') if stage == 'first' else ('F', 'R', 'T')
                        for cls in classes:
                            original = base if cls in ('D', 'F') else base[base.transition == cls]
                            view = kept if cls in ('D', 'F') else kept[kept.transition == cls]
                            meta = dict(mask=mask, organism=organism, scheme=scheme, stage=stage, threshold=threshold, transition=cls)
                            p0 = original[original.species == organism]
                            p1 = view[view.species == organism]
                            a0 = original[original.species.isin(('H', 'Y', 'E')) & (original.species != organism)]
                            a1 = view[view.species.isin(('H', 'Y', 'E')) & (view.species != organism)]
                            coverage.append(dict(**meta, accepted_baseline=len(original), accepted_retained=len(view), positive_raw_baseline=int(positive(original['Precursor.Quantity']).sum()), positive_raw_retained=int(positive(view['Precursor.Quantity']).sum()), precursor_rows_baseline=original.key.nunique(), precursor_rows_retained=view.key.nunique(), selected_filter_rows=len(selected), designed_runs=6, possible_selected_entries=len(selected) * 6, present_events_baseline=len(p0), present_events_retained=len(p1), present_positive_baseline=int(positive(p0['Precursor.Quantity']).sum()), present_positive_retained=int(positive(p1['Precursor.Quantity']).sum()), present_rows_baseline=p0.key.nunique(), present_rows_retained=p1.key.nunique(), absent_events_baseline=len(a0), absent_events_retained=len(a1), absent_positive_baseline=int(positive(a0['Precursor.Quantity']).sum()), absent_positive_retained=int(positive(a1['Precursor.Quantity']).sum()), absent_rows_baseline=a0.key.nunique(), absent_rows_retained=a1.key.nunique()))
                            for scope, part in [('group', view)] + [(r, view[view['sample'] == r]) for r in axis]:
                                for unit in ('events', 'rows'):
                                    unitdata = part if unit == 'events' else part.drop_duplicates('key')
                                    A = int((unitdata.species.isin(('H', 'Y', 'E')) & (unitdata.species != organism)).sum())
                                    P = int((unitdata.species == organism).sum())
                                    rid = len(records)
                                    records.append(dict(**meta, scope=scope, unit=unit, A=A, P=P, n=A + P, false_fraction=A / (A + P) if A + P else np.nan, status='estimable' if A + P else 'not_estimable', shared=int((unitdata.species == 'shared').sum()), shared_IL=int((unitdata.species == 'shared_IL').sum()), unknown=int((unitdata.species == 'unknown').sum()), accepted_count=len(unitdata)))
                                    exclusive = part[part.species.isin(('H', 'Y', 'E'))]
                                    for label, group in exclusive.groupby('species', sort=True):
                                        side = 0 if label != organism else 1
                                        if unit == 'events':
                                            for sample, k in group[['sample', 'key']].itertuples(index=False, name=None):
                                                contributions.append((rid, side, keyidx[k], 1 << runidx[sample], 0))
                                        else:
                                            for k, kg in group.groupby('key', sort=True):
                                                bit = sum((1 << runidx[r] for r in sorted(set(kg['sample']))))
                                                contributions.append((rid, side, keyidx[k], bit, 1))
    table = pd.DataFrame(records)
    cov = pd.DataFrame(coverage)
    if len(cov):
        for name, a, b in [('present_event_coverage', 'present_events_retained', 'present_events_baseline'), ('present_positive_coverage', 'present_positive_retained', 'present_positive_baseline'), ('present_row_coverage', 'present_rows_retained', 'present_rows_baseline'), ('absent_event_retention', 'absent_events_retained', 'absent_events_baseline'), ('absent_positive_retention', 'absent_positive_retained', 'absent_positive_baseline'), ('absent_row_retention', 'absent_rows_retained', 'absent_rows_baseline')]:
            cov[name] = np.divide(cov[a], cov[b], out=np.full(len(cov), np.nan), where=cov[b] > 0)
    pooled_indices = []
    if len(table):
        groups = ['mask', 'scheme', 'stage', 'threshold', 'transition', 'unit']
        for values, g in table[table.scope == 'group'].groupby(groups, sort=True):
            if len(g) != 3:
                continue
            row = dict(zip(groups, values))
            row.update(organism='ALL', scope='pooled_organism_key' if row['unit'] == 'rows' else 'pooled_events')
            for c in ['A', 'P', 'n', 'shared', 'shared_IL', 'unknown', 'accepted_count']:
                row[c] = int(g[c].sum())
            row['false_fraction'] = row['A'] / row['n'] if row['n'] else np.nan
            row['status'] = 'estimable' if row['n'] else 'not_estimable'
            pooled_indices.append((len(table), list(g.index)))
            table.loc[len(table)] = row
    if draws and len(table):
        rng = np.random.default_rng(SEED)
        out = {b: np.full((draws, len(table)), np.nan) for b in ('recipient', 'precursor', 'two_way')}
        if contributions:
            c = np.array(contributions, dtype=np.int64)
            bits, binv = np.unique(c[:, 3], return_inverse=True)
            inc = (bits[:, None] >> np.arange(18) & 1).astype(np.int8)
            reduce = sparse.csr_matrix((np.ones(len(c)), (c[:, 0] * 2 + c[:, 1], np.arange(len(c)))), shape=(2 * len(table), len(c)))
            strata = [np.array([runidx[r] for r in singles[singles.present_species == sp]['sample']], int) for sp in ('H', 'Y', 'E')]
            for start in range(0, draws, 16):
                n = min(16, draws - start)
                kw = rng.multinomial(len(keylist), np.full(len(keylist), 1 / len(keylist)), size=n)
                rw = np.zeros((n, 18), int)
                for ids in strata:
                    rw[:, ids] = rng.multinomial(6, np.full(6, 1 / 6), size=n)
                rsum = rw @ inc.T
                rweights = np.where(c[:, 4][None, :] == 1, rsum[:, binv] > 0, rsum[:, binv])
                for model in out:
                    weight = (rweights if model != 'precursor' else np.ones_like(rweights)) * (kw[:, c[:, 2]] if model != 'recipient' else 1)
                    ap = (reduce @ weight.T).T.reshape(n, len(table), 2)
                    for idx, children in pooled_indices:
                        ap[:, idx, :] = ap[:, children, :].sum(axis=1)
                    den = ap.sum(axis=2)
                    out[model][start:start + n] = np.divide(ap[:, :, 0], den, out=np.full(den.shape, np.nan), where=den > 0)
        for model, values in out.items():
            ci = np.array([interval(values[:, i]) for i in range(len(table))])
            for idx, label in enumerate(('ci_low', 'ci_high', 'undefined_draws', 'infinite_draws')):
                table[model + '_' + label] = ci[:, idx]
    return (table, cov, pd.DataFrame(membership, columns=['mask', 'organism', 'scheme', 'stage', 'selector_stage', 'threshold', 'key', 'observed_runs', 'designed_runs', 'selected']), pd.DataFrame(unavailable, columns=['mask', 'organism', 'scheme', 'stage', 'status', 'reason']))

def fc(e, s, status, keys, draws=2000, *, guard=None, mixture_date='20250116'):
    producer = Path(__file__).with_name('fold_change_replication.py') if mixture_date == '20250109' else __file__
    require_guard(guard, producer=producer, frames=(e, s, status))
    axis = s[(s.kind == 'mixture') & s['group'].isin(['5H4Y1E', '5H1Y4E'])]
    lookup = axis.set_index('sample')['group'].to_dict()
    subset = e[e['sample'].isin(lookup) & e.species.isin(['Y', 'E'])].copy()
    records = []
    unavailable = []
    for mask in MASKS:
        ok = status[status['sample'].isin(lookup) & (status['mask'] == mask)].status.eq('complete').all()
        if not ok:
            unavailable.append(dict(mask=mask, status='not_estimable', reason='required extreme-mixture mask/stage unavailable'))
            continue
        m = subset[subset['mask'] == mask]
        for sp in ('Y', 'E'):
            universe = sorted(m[m.species == sp].key.unique())
            if not universe:
                unavailable.append(dict(mask=mask, species=sp, status='not_estimable', reason='no accepted species-exclusive precursor keys in extreme mixtures'))
            for context in ('D', 'F', 'R', 'T'):
                part = m[(m.stage == ('first' if context == 'D' else 'final')) & (m.species == sp)]
                if context in ('R', 'T'):
                    part = part[part.transition == context]
                for quantity in ('Precursor.Normalised', 'Precursor.Quantity'):
                    p = part[positive(part[quantity])].copy()
                    p['group'] = p['sample'].map(lookup)
                    p['log2'] = np.log2(p[quantity])
                    agg = p.groupby(['key', 'group']).log2.agg(['count', 'median'])
                    for k in universe:
                        hi = agg.loc[k, '5H4Y1E'] if (k, '5H4Y1E') in agg.index else None
                        lo = agg.loc[k, '5H1Y4E'] if (k, '5H1Y4E') in agg.index else None
                        nh = 0 if hi is None else int(hi['count'])
                        nl = 0 if lo is None else int(lo['count'])
                        eligible = nh >= 2 and nl >= 2
                        observed = float(hi['median'] - lo['median']) if eligible else np.nan
                        expected = 2.0 if sp == 'Y' else -2.0
                        records.append(dict(mask=mask, species=sp, context=context, quantity=quantity, key=k, batch=mixture_date, n_high=nh, n_low=nl, eligible=eligible, observed_log2fc=observed, expected_log2fc=expected, error=observed - expected, coverage_population='own_context'))
    rowcols = ['mask', 'species', 'context', 'quantity', 'key', 'batch', 'n_high', 'n_low', 'eligible', 'observed_log2fc', 'expected_log2fc', 'error', 'coverage_population']
    rows = pd.DataFrame(records, columns=rowcols)
    summaries = []
    if not len(rows):
        return (rows, pd.DataFrame(columns=['mask', 'species', 'quantity', 'context', 'population', 'eligible_keys', 'mean_bias', 'median_bias', 'MAE', 'RMSE']), pd.DataFrame(unavailable, columns=['mask', 'species', 'status', 'reason']))
    for (mask, sp, quantity), g in rows.groupby(['mask', 'species', 'quantity'], sort=True):
        pivot = g.pivot(index='key', columns='context', values='error').reindex(columns=['D', 'F', 'R', 'T'])
        plans = [(c, 'own_context', pivot[c].notna()) for c in ['D', 'F', 'R', 'T']]
        for a, b in [('D', 'F'), ('R', 'T')]:
            shared = pivot[a].notna() & pivot[b].notna()
            plans.extend([(a, a + '_' + b + '_intersection', shared), (b, a + '_' + b + '_intersection', shared)])
        rng = np.random.default_rng(SEED)
        for context, pop, selection in plans:
            values = pivot.loc[selection, context].to_numpy(float)
            expected = 2.0 if sp == 'Y' else -2.0
            record = dict(mask=mask, species=sp, quantity=quantity, context=context, population=pop, eligible_keys=len(values), universe_keys=len(pivot), unavailable_keys=len(pivot) - len(values), mean_bias=float(values.mean()) if len(values) else np.nan, median_bias=float(np.median(values)) if len(values) else np.nan, MAE=float(np.abs(values).mean()) if len(values) else np.nan, RMSE=float(np.sqrt(np.square(values).mean())) if len(values) else np.nan, expected_log2fc=expected)
            record['mean_observed_log2fc'] = record['mean_bias'] + expected
            record['median_observed_log2fc'] = record['median_bias'] + expected
            if draws:
                rng = np.random.default_rng(SEED)
                allvalues = pivot[context].where(selection).to_numpy(float)
                boot = {k: [] for k in ['mean_bias', 'median_bias', 'MAE', 'RMSE']}
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    for start in range(0, draws, 16):
                        n = min(16, draws - start)
                        selected = allvalues[rng.integers(0, len(pivot), size=(n, len(pivot)))]
                        boot['mean_bias'].extend(np.nanmean(selected, axis=1))
                        boot['median_bias'].extend(np.nanmedian(selected, axis=1))
                        boot['MAE'].extend(np.nanmean(np.abs(selected), axis=1))
                        boot['RMSE'].extend(np.sqrt(np.nanmean(selected ** 2, axis=1)))
                for metric, b in boot.items():
                    for tag, v in zip(['ci_low', 'ci_high', 'undefined_draws', 'infinite_draws'], interval(b)):
                        record[metric + '_' + tag] = v
            summaries.append(record)
    return (rows, pd.DataFrame(summaries), pd.DataFrame(unavailable, columns=['mask', 'species', 'status', 'reason']))

def output_frame(path, frame):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, lineterminator='\n', na_rep='NA')
