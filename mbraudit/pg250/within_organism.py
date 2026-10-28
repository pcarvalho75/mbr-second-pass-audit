"""Within-organism co-searches: six single-organism runs analysed alone (carried from validation_v5/s1_adapter.py)."""
from collections import Counter
from pathlib import Path
import hashlib, json, math
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix


HERE = Path(__file__).resolve().parent


KEY = ['Modified.Sequence', 'Precursor.Charge']


EVENT = ['sample'] + KEY


MASKS = {'F0': ['Q.Value'], 'F3': ['Q.Value', 'Global.Q.Value', 'Lib.Q.Value'], 'F4': ['Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value', 'Lib.PG.Q.Value']}


SPECIES = ['H', 'Y', 'E', 'shared', 'shared_IL', 'unknown', 'quarantined']


BITS_SPECIES = ['H', 'Y', 'E']


TERCILES = ['low', 'middle', 'high', 'missing_quantity']


BINS = ['zero', 'positive_0_20', 'positive_20_40', 'positive_40_60', 'positive_60_100']


SEED = 20261018


DRAWS = 2000


BATCH = 16


INPUT_ROLES = {'first', 'final', 'empirical', 'predicted_species', 'design', 'exports'}


def table(path, frame):
    return _base().table(path, frame)


def fraction(a, p):
    return _base().fraction(a, p)


def ratio(a, p):
    return _base().ratio(a, p)


def interval(values):
    return _base().interval(values)


def support_bin(k, denominator=5):
    return _base().support_bin(k, denominator)


def quantity_terciles(values):
    return _base().quantity_terciles(values)


def canonical_sample(value):
    return _base().canonical_sample(value)


def scope(organism):
    return {'analysis': 'within_organism', 'organism': organism, 'cohort_size': 6, 'support_denominator': 5}


def design_check(design, organism, synthetic=False):
    required = ['sample', 'kind', 'group', 'filename_date', 'present_species']
    if any((k not in design for k in required)):
        raise ValueError('Incomplete within-organism design')
    if len(design) != 6 or design['sample'].duplicated().any() or design['sample'].isna().any():
        raise ValueError('Exactly six unique designed axes required')
    if not design.kind.eq('single').all() or not design.present_species.eq(organism).all() or (not design.group.eq(organism).all()):
        raise ValueError('One nominal-organism stratum of six single runs required')
    if synthetic and (not design['sample'].str.startswith('SYN_').all()):
        raise ValueError('Synthetic sample identities required')
    return design


def interface_check(events, status, design):
    expected = {(s, m, t) for s in design['sample'] for m in MASKS for t in ['first', 'final']}
    if status.duplicated(['sample', 'mask', 'stage']).any() or set(map(tuple, status[['sample', 'mask', 'stage']].to_numpy())) != expected:
        raise ValueError('Every designed run/mask/stage status must be explicit')
    if not status.status.isin(['complete', 'unavailable']).all():
        raise ValueError('Invalid stage status')
    if events.duplicated(['mask', 'stage'] + EVENT).any() or not set(events['sample']) <= set(design['sample']):
        raise ValueError('Duplicate or unknown accepted event')
    if not set(events.species) <= set(SPECIES) or not set(events['mask']) <= set(MASKS) or (not set(events.stage) <= {'first', 'final'}):
        raise ValueError('Unknown accepted-event category')


def memberships(events, status, design, empirical, organism):
    design_check(design, organism)
    interface_check(events, status, design)
    members, available = _memberships_frozen_delta(events, status, design, empirical)
    if members.empty:
        for name in ['Precursor.Quantity', 'Precursor.Normalised']:
            members[name] = pd.Series(index=members.index, dtype=float)
    return (members, available)


def endpoint_summary(members, design, available, organism):
    design_check(design, organism)
    summary, a, p, pairs, keys = _endpoint_summary_frozen_delta(members, design, available)
    if summary.empty:
        columns = ['endpoint_id', 'endpoint', 'mask', 'class', 'tercile', 'support_bin', 'scope', 'group', 'events', 'A', 'P', 'n', 'absent_fraction', 'distinct_precursors', 'distinct_stripped_peptides', 'distinct_accession_groups', 'positive_quantity_events']
        columns += ['excluded_' + sp for sp in SPECIES if sp not in BITS_SPECIES]
        summary = pd.DataFrame(columns=columns)
    return (summary, a, p, pairs, keys)


def bootstrap(summary, a, p, pairs, keys, design, out, draws=2000):
    design_check(design, design.present_species.iloc[0])
    return _bootstrap_frozen_delta(summary, a, p, pairs, keys, design, out, draws)


def _memberships_frozen_delta(events, status, design, empirical):
    singles = design[design.kind.eq('single')]
    single_names = set(singles['sample'])
    result = []
    availability = []
    for mask in MASKS:
        axes = status[status['mask'].eq(mask)]
        pair_ok = axes[axes['sample'].isin(single_names)].status.eq('complete').all()
        support_ok = axes[axes.stage.eq('first')].status.eq('complete').all()
        availability.append({'mask': mask, 'T1_available': bool(pair_ok), 'T2_available': bool(pair_ok and support_ok), 'support_denominator': 5, 'loss_library_split_available': empirical is not None})
        if not pair_ok:
            continue
        frame = events[events['mask'].eq(mask)]
        first = frame[frame.stage.eq('first')]
        final = frame[frame.stage.eq('final')]
        support = first.groupby(KEY).size().to_dict()
        for sample in singles['sample']:
            d = first[first['sample'].eq(sample)].set_index(KEY, drop=False)
            f = final[final['sample'].eq(sample)].set_index(KEY, drop=False)
            D = set(d.index)
            F = set(f.index)
            R = D & F
            T = F - D
            L = D - F
            assert len(D) == len(R) + len(L) and len(F) == len(R) + len(T)
            for label, keys, source in [('D', D, d), ('F', F, f), ('R', R, f), ('T', T, f), ('L', L, d)]:
                if not keys:
                    continue
                chosen = source.loc[sorted(keys)].copy()
                chosen['class'] = label
                if label == 'L':
                    chosen['library_state'] = ['unavailable' if empirical is None else 'in_library' if key in empirical else 'outside_library' for key in chosen.index]
                else:
                    chosen['library_state'] = 'not_applicable'
                chosen['tercile'] = 'not_applicable'
                chosen['support_bin'] = 'not_applicable'
                chosen['support_k'] = pd.NA
                if label == 'T':
                    chosen['tercile'] = quantity_terciles(chosen['Precursor.Quantity'])
                    if support_ok:
                        chosen['support_k'] = [support.get(key, 0) - int(key in D) for key in chosen.index]
                        chosen['support_bin'] = [support_bin(int(k), 5) for k in chosen.support_k]
                    else:
                        chosen['support_bin'] = 'unavailable'
                result.append(chosen.reset_index(drop=True))
    columns = list(events.columns) + ['class', 'library_state', 'tercile', 'support_bin', 'support_k']
    return (pd.concat(result, ignore_index=True) if result else pd.DataFrame(columns=columns), pd.DataFrame(availability))


def _endpoint_summary_frozen_delta(members, design, available):
    """Build every empty axis explicitly and a shared sparse resampling membership."""
    singles = design[design.kind.eq('single')]
    present = singles.set_index('sample').present_species.to_dict()
    names = singles['sample'].tolist()
    keys = sorted(set(map(tuple, members[KEY].to_numpy())))
    pairkeys = sorted(set(map(tuple, members[EVENT].to_numpy())), key=lambda x: (names.index(x[0]), x[1], x[2]))
    pairindex = {key: i for i, key in enumerate(pairkeys)}
    keyindex = {key: i for i, key in enumerate(keys)}
    pairs = pd.DataFrame(pairkeys, columns=EVENT)
    pairs['run_index'] = [names.index(x) for x in pairs['sample']]
    pairs['key_index'] = [keyindex[tuple(x)] for x in pairs[KEY].to_numpy()]
    scopes = [('run', sample, [sample]) for sample in names] + [('organism', singles.present_species.iloc[0], names)]
    records = []
    ar = []
    ac = []
    pr = []
    pc = []
    contribution = []
    for mask in MASKS:
        av = available[available['mask'].eq(mask)].iloc[0]
        for endpoint in ['T1', 'T2']:
            if not av[endpoint + '_available']:
                continue
            cells = [(cls, 'all', 'all') for cls in ['D', 'F', 'R', 'T', 'L', 'L_in_library', 'L_outside_library'] if not cls.startswith('L_') or av.loss_library_split_available] if endpoint == 'T1' else [('T', q, b) for q in TERCILES for b in BINS]
            base = members[members['mask'].eq(mask)]
            for cls, tercile, bin_name in cells:
                selected = base[base['class'].eq(cls if not cls.startswith('L_') else 'L')]
                if cls.startswith('L_'):
                    selected = selected[selected.library_state.eq(cls[2:])]
                if endpoint == 'T2':
                    selected = selected[selected.tercile.eq(tercile) & selected.support_bin.eq(bin_name)]
                for scope, group, samples in scopes:
                    data = selected[selected['sample'].isin(samples)]
                    i = len(records)
                    counts = Counter(data.species)
                    A = P = 0
                    for event in data.itertuples(index=False, name=None):
                        row = dict(zip(data.columns, event))
                        sp = row['species']
                        if sp in BITS_SPECIES:
                            pair = (row['sample'], row[KEY[0]], row[KEY[1]])
                            if sp == present[row['sample']]:
                                P += 1
                                pr.append(i)
                                pc.append(pairindex[pair])
                            else:
                                A += 1
                                ar.append(i)
                                ac.append(pairindex[pair])
                    records.append({'endpoint_id': i, 'endpoint': endpoint, 'mask': mask, 'class': cls, 'tercile': tercile, 'support_bin': bin_name, 'scope': scope, 'group': group, 'events': len(data), 'A': A, 'P': P, 'n': A + P, 'absent_fraction': float(fraction(A, P)), 'distinct_precursors': len(data[KEY].drop_duplicates()), 'distinct_stripped_peptides': data['Stripped.Sequence'].replace('', pd.NA).nunique(), 'distinct_accession_groups': data['Protein.Ids'].replace('', pd.NA).nunique(), 'positive_quantity_events': int((data['Precursor.Quantity'].gt(0) & np.isfinite(data['Precursor.Quantity'])).sum()), **{'excluded_' + s: int(counts.get(s, 0)) for s in SPECIES if s not in BITS_SPECIES}})
    shape = (len(records), len(pairkeys))
    a = coo_matrix((np.ones(len(ar), dtype=np.int64), (ar, ac)), shape=shape).tocsr()
    p = coo_matrix((np.ones(len(pr), dtype=np.int64), (pr, pc)), shape=shape).tocsr()
    return (pd.DataFrame(records), a, p, pairs, keys)


def _bootstrap_frozen_delta(summary, a, p, pairs, keys, design, out, draws=2000):
    """One shared recipient/key draw per replicate across all filters and endpoints."""
    if draws != 2000:
        raise ValueError('Frozen draw count is 2000')
    singles = design[design.kind.eq('single')]
    names = singles['sample'].tolist()
    strata = [np.arange(6, dtype=int)]
    ri = pairs.run_index.to_numpy(dtype=int)
    ki = pairs.key_index.to_numpy(dtype=int)
    rng = np.random.Generator(np.random.PCG64(SEED))
    nkeys = len(keys)
    ledger = {scheme: np.zeros((2, len(summary), draws), dtype=np.int64) for scheme in ['recipient', 'precursor', 'two_way']}
    for start in range(0, draws, 16):
        width = min(16, draws - start)
        rw = np.zeros((6, width), dtype=np.int64)
        kw = np.zeros((nkeys, width), dtype=np.int64)
        for column in range(width):
            for group in strata:
                rw[:, column] += np.bincount(rng.choice(group, len(group), replace=True), minlength=6)
            if nkeys:
                kw[:, column] = np.bincount(rng.integers(nkeys, size=nkeys), minlength=nkeys)
        weights = {'recipient': rw[ri], 'precursor': kw[ki], 'two_way': rw[ri] * kw[ki]}
        for scheme, weight in weights.items():
            ledger[scheme][0, :, start:start + width] = a @ weight
            ledger[scheme][1, :, start:start + width] = p @ weight
    np.savez_compressed(out / 'bootstrap_counts.npz', **ledger)
    np.savez_compressed(out / 'bootstrap_fractions.npz', **{scheme: fraction(counts[0], counts[1]) for scheme, counts in ledger.items()})
    intervals = []
    for scheme, counts in ledger.items():
        fractions = fraction(counts[0], counts[1])
        for index, values in enumerate(fractions):
            low, high, bad, inf, state = interval(values)
            intervals.append({'endpoint_id': index, 'scheme': scheme, 'draws': draws, 'point': summary.iloc[index].absent_fraction, 'low95': low, 'high95': high, 'undefined_draws': bad, 'infinite_draws': inf, 'interval_status': state})
    table(out / 'fraction_intervals.tsv', pd.DataFrame(intervals))
    return ledger


def planned_cells(design, available, summary):
    organism = design.present_species.iloc[0]
    scopes = [('run', s) for s in design['sample']] + [('organism', organism)]
    seen = {(r.endpoint, r['mask'], r['class'], r.tercile, r.support_bin, r.scope, r.group): int(r.endpoint_id) for _, r in summary.iterrows()}
    rows = []
    for mask in MASKS:
        av = available[available['mask'].eq(mask)].iloc[0]
        for endpoint in ['T1', 'T2']:
            cells = [(c, 'all', 'all') for c in ['D', 'F', 'R', 'T', 'L', 'L_in_library', 'L_outside_library']] if endpoint == 'T1' else [('T', q, b) for q in TERCILES for b in BINS]
            for cls, tercile, bin_name in cells:
                ok = bool(av[endpoint + '_available']) and (not cls.startswith('L_') or bool(av.loss_library_split_available))
                for scope_name, group in scopes:
                    key = (endpoint, mask, cls, tercile, bin_name, scope_name, group)
                    rows.append({'endpoint': endpoint, 'mask': mask, 'class': cls, 'tercile': tercile, 'support_bin': bin_name, 'scope': scope_name, 'group': group, 'available': ok, 'endpoint_id': seen.get(key, pd.NA), 'status': 'available' if ok else 'unavailable_required_stage_q_field_or_empirical', 'display': 'primary' if mask == 'F3' and (endpoint == 'T1' or tercile == 'high') else 'sensitivity'})
    return pd.DataFrame(rows)


def interpretation_labels(out):
    path = out / 'contrast_intervals.tsv'
    if not path.is_file():
        return
    frame = pd.read_csv(path, sep='\t')
    finite = np.isfinite(frame.low95) & np.isfinite(frame.high95)
    frame['interpretation_available'] = finite
    frame['interpretation_status'] = np.where(finite, 'finite_population_conditional_interval', 'unavailable_or_infinite_interval')
    primary = frame['mask'].eq('F3') & frame['scheme'].eq('two_way') & (frame.endpoint.eq('T1') | frame.tercile.eq('high'))
    frame['display'] = np.where(primary, 'primary_interval', 'sensitivity_interval')
    table(out / 'contrast_interpretation.tsv', frame)


def _base():
    from mbraudit.pg250 import endpoints as exact_negative_endpoints
    return exact_negative_endpoints
