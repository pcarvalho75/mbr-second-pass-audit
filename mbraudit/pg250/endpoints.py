"""Exact-negative accounting of the 250 pg dataset: export, classes, A/P counts, bootstrap and contrasts (carried from validation_v5/exact_negative_endpoints.py)."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from mbraudit.library.opportunity_audit import strip_mods
KEY = ['Modified.Sequence', 'Precursor.Charge']
EVENT = ['sample'] + KEY
MASKS = {'F0': ['Q.Value'], 'F3': ['Q.Value', 'Global.Q.Value', 'Lib.Q.Value'], 'F4': ['Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value', 'Lib.PG.Q.Value']}
SPECIES = ['H', 'Y', 'E', 'shared', 'shared_IL', 'unknown', 'quarantined']
TERCILES = ['low', 'middle', 'high', 'missing_quantity']
BINS = ['zero', 'positive_0_20', 'positive_20_40', 'positive_40_60', 'positive_60_100']
SEED = 20261018

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def table(path, frame):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, na_rep='NA', lineterminator='\n')

def json_write(path, value):
    Path(path).write_bytes((json.dumps(value, indent=2, allow_nan=False) + '\n').encode())

def design_check(design, synthetic=False, *, mixture_date='20250116'):
    required = ['sample', 'kind', 'group', 'filename_date', 'present_species']
    if any((x not in design for x in required)):
        raise ValueError('Incomplete sample metadata')
    if len(design) != 54 or design['sample'].duplicated().any():
        raise ValueError('Expected 54 unique designed axes')
    if synthetic and (not design['sample'].str.startswith('SYN_').all()):
        raise ValueError('Synthetic sample names required')
    singles = design[design.kind.eq('single')]
    if len(singles) != 18 or singles.present_species.value_counts().to_dict() != {'H': 6, 'Y': 6, 'E': 6}:
        raise ValueError('Expected six single-organism recipients per species')
    mix = design[design.kind.eq('mixture')]
    if len(mix) != 30 or set(mix.filename_date.astype(str)) != {str(mixture_date)}:
        raise ValueError('Expected 30 mixtures in the fixed batch')
    if mix.group.value_counts().to_dict() != {x: 6 for x in ['5H1Y4E', '5H2Y3E', '4H2Y2E', '5H3Y2E', '5H4Y1E']}:
        raise ValueError('Mixture groups disagree')
    if int(design.kind.eq('blank').sum()) != 6:
        raise ValueError('Expected six blank axes')
    return singles

def canonical_sample(value):
    value = str(value).replace('\\', '/').rsplit('/', 1)[-1]
    for suffix in ['.d.dia', '.raw.dia', '.mzML', '.raw', '.dia', '.d']:
        if value.lower().endswith(suffix.lower()):
            return value[:-len(suffix)]
    return value

def support_bin(k, denominator=53):
    if not 0 <= k <= denominator:
        raise ValueError('Invalid support count')
    if k == 0:
        return 'zero'
    for multiple, label in [(1, BINS[1]), (2, BINS[2]), (3, BINS[3])]:
        if 5 * k <= multiple * denominator:
            return label
    return BINS[4]

def quantity_terciles(values):
    values = pd.to_numeric(values, errors='coerce')
    good = values.gt(0) & np.isfinite(values)
    result = pd.Series('missing_quantity', index=values.index, dtype='str')
    m = int(good.sum())
    if m:
        u = (values[good].rank(method='average') - 0.5) / m
        result.loc[good] = np.where(u <= 1 / 3, 'low', np.where(u <= 2 / 3, 'middle', 'high'))
    return result

def duplicate_audit(frame, stage, out, *, stop_on_conflict=True):
    duplicated = frame.duplicated(EVENT, keep=False)
    ledger = frame.loc[duplicated].copy()
    table(out / (stage + '_duplicate_rows.tsv'), ledger)
    unique = frame.drop_duplicates()
    if stop_on_conflict and unique.duplicated(EVENT).any():
        raise ValueError('Conflicting duplicate run/key rows; see duplicate audit')
    return (unique, {'stage': stage, 'duplicate_rows': int(duplicated.sum()), 'identical_rows_removed': len(frame) - len(unique)})

def adapt_reports(config, design, out):
    """Produce the accepted-event/run-status interface, retaining zero-observation axes."""
    inputs = config['inputs']
    annotation = pd.read_parquet(inputs['predicted_species'])
    if annotation.duplicated(KEY).any():
        raise ValueError('Predicted map must have unique frozen keys')
    if not set(annotation.species).issubset(SPECIES):
        raise ValueError('Unknown species label')
    annotation = annotation[KEY + ['species', 'Stripped.Sequence', 'Protein.Ids']].rename(columns={'Protein.Ids': 'frozen_Protein.Ids', 'Stripped.Sequence': 'frozen_Stripped.Sequence'})
    exports = pd.read_csv(inputs['exports'], sep='\t')
    expected = {(sample, stage) for sample in design['sample'] for stage in ['first', 'final']}
    if exports.duplicated(['sample', 'stage']).any() or set(map(tuple, exports[['sample', 'stage']].to_numpy())) != expected:
        raise ValueError('Missing/duplicate designed export receipts')
    if not exports.status.isin(['complete', 'unavailable']).all():
        raise ValueError('Invalid export state')
    accepted = []
    statuses = []
    q_audit = []
    duplicates = []
    decoys = []
    decoy_duplicates = []
    for stage in ['first', 'final']:
        raw = pd.read_parquet(inputs[stage])
        if 'sample' not in raw:
            field = 'Run' if 'Run' in raw else 'File.Name'
            if field not in raw:
                raise ValueError('Missing report run identity')
            raw['sample'] = raw[field].map(canonical_sample)
        if not set(raw['sample']).issubset(set(design['sample'])):
            raise ValueError('Unexpected report runs')
        if any((x not in raw for x in KEY + ['Decoy'])):
            raise ValueError('Missing identity/decoy fields')
        charges = pd.to_numeric(raw['Precursor.Charge'], errors='raise')
        if not np.isfinite(charges).all() or not (charges == np.floor(charges)).all():
            raise ValueError('Invalid integer charge')
        raw['Precursor.Charge'] = charges.astype(int)
        decoy_rows = raw[raw.Decoy.ne(0)].copy()
        raw = raw[raw.Decoy.eq(0)].copy()
        table(out / (stage + '_decoy_rows.tsv'), decoy_rows)
        decoy_raw_counts = decoy_rows['sample'].value_counts()
        decoy_rows, decoy_dup = duplicate_audit(decoy_rows, stage + '_decoy', out, stop_on_conflict=False)
        conflicts = decoy_rows[decoy_rows.duplicated(EVENT, keep=False)]
        decoy_duplicates.append({**decoy_dup, 'conflicting_groups': len(conflicts[EVENT].drop_duplicates())})
        raw, dup = duplicate_audit(raw, stage, out)
        duplicates.append(dup)
        q_source = pd.concat([raw, decoy_rows]).sort_index(kind='stable')
        raw = raw.merge(annotation, on=KEY, how='left', validate='many_to_one')
        raw['species'] = raw.species.fillna('unknown')
        for col in ['Precursor.Quantity', 'Precursor.Normalised']:
            if col not in raw:
                raw[col] = np.nan
            raw[col] = pd.to_numeric(raw[col], errors='coerce')
        for field in MASKS['F4']:
            if field in raw:
                q = pd.to_numeric(q_source[field], errors='coerce')
                raw[field] = pd.to_numeric(raw[field], errors='coerce')
                decoy_rows[field] = pd.to_numeric(decoy_rows[field], errors='coerce')
                finite = q[np.isfinite(q)]
                q_audit.append({'stage': stage, 'field': field, 'present': True, 'rows': len(q), 'nonfinite': int((~np.isfinite(q)).sum()), 'below_zero': int(q.lt(0).sum()), 'above_one': int(q.gt(1).sum()), 'unique_finite': finite.nunique(), 'all_zero': bool(len(finite) == len(q) and len(q) > 0 and finite.eq(0).all()), 'constant_finite': bool(len(finite) > 0 and finite.nunique() == 1)})
            else:
                q_audit.append({'stage': stage, 'field': field, 'present': False, 'rows': len(q_source)})
        for mask, fields in MASKS.items():
            available = all((x in raw for x in fields))
            keep = pd.Series(available, index=raw.index)
            decoy_keep = pd.Series(available, index=decoy_rows.index)
            if available:
                for field in fields:
                    keep &= raw[field].le(0.01) & np.isfinite(raw[field])
                    decoy_keep &= decoy_rows[field].le(0.01) & np.isfinite(decoy_rows[field])
            axes = exports[exports.stage.eq(stage)].set_index('sample').status
            for sample in design['sample']:
                status = 'complete' if available and axes[sample] == 'complete' else 'unavailable'
                statuses.append({'sample': sample, 'stage': stage, 'mask': mask, 'status': status, 'reason': '' if status == 'complete' else 'missing_q_field' if not available else 'export_unavailable'})
            invalid = set(axes[axes.ne('complete')].index)
            if raw['sample'].isin(invalid).any() or decoy_rows['sample'].isin(invalid).any():
                raise ValueError('Rows present for an unavailable export')
            dec = decoy_rows[decoy_keep]
            for sample in design['sample']:
                decoys.append({'sample': sample, 'stage': stage, 'mask': mask, 'accepted_decoy_rows': int(dec['sample'].eq(sample).sum()), 'raw_decoy_rows': int(decoy_raw_counts.get(sample, 0))})
            selected = raw[keep].copy()
            selected['mask'] = mask
            selected['stage'] = stage
            selected['Protein.Ids'] = selected['frozen_Protein.Ids'].fillna('')
            selected['Stripped.Sequence'] = [frozen if pd.notna(frozen) else strip_mods(modified) for frozen, modified in zip(selected['frozen_Stripped.Sequence'], selected['Modified.Sequence'])]
            accepted.append(selected[['sample'] + KEY + ['species', 'mask', 'stage', 'Precursor.Quantity', 'Precursor.Normalised', 'Stripped.Sequence', 'Protein.Ids']])
    events = pd.concat(accepted, ignore_index=True)
    status = pd.DataFrame(statuses)
    if events.duplicated(['mask', 'stage'] + EVENT).any():
        raise ValueError('Duplicate accepted event')
    exported = events.copy()
    exported['annotation_status'] = np.where(exported.species.eq('quarantined'), 'quarantined', 'frozen_or_unknown')
    exported.loc[exported.species.eq('quarantined'), 'species'] = 'unknown'
    for name, frame in [('accepted_events', exported), ('run_status', status), ('q_field_audit', pd.DataFrame(q_audit)), ('duplicate_audit', pd.DataFrame(duplicates)), ('decoy_audit', pd.DataFrame(decoys)), ('decoy_duplicate_audit', pd.DataFrame(decoy_duplicates))]:
        table(out / (name + '.tsv'), frame)
    empirical = None
    if inputs.get('empirical'):
        empirical = pd.read_parquet(inputs['empirical'])
        if 'Decoy' not in empirical:
            raise ValueError('Missing empirical-library decoy field')
        empirical = empirical[empirical.Decoy.eq(0)]
        empirical = set(map(tuple, empirical[KEY].drop_duplicates().to_numpy()))
    return (events, status, empirical)

def memberships(events, status, design, empirical):
    singles = design[design.kind.eq('single')]
    single_names = set(singles['sample'])
    result = []
    availability = []
    for mask in MASKS:
        axes = status[status['mask'].eq(mask)]
        pair_ok = axes[axes['sample'].isin(single_names)].status.eq('complete').all()
        support_ok = axes[axes.stage.eq('first')].status.eq('complete').all()
        availability.append({'mask': mask, 'T1_available': bool(pair_ok), 'T2_available': bool(pair_ok and support_ok), 'support_denominator': 53, 'loss_library_split_available': empirical is not None})
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
                        chosen['support_bin'] = [support_bin(int(k)) for k in chosen.support_k]
                    else:
                        chosen['support_bin'] = 'unavailable'
                result.append(chosen.reset_index(drop=True))
    columns = list(events.columns) + ['class', 'library_state', 'tercile', 'support_bin', 'support_k']
    return (pd.concat(result, ignore_index=True) if result else pd.DataFrame(columns=columns), pd.DataFrame(availability))

def fraction(a, p):
    a = np.asarray(a, dtype=float)
    n = a + np.asarray(p, dtype=float)
    return np.divide(a, n, out=np.full(np.broadcast(a, n).shape, np.nan), where=n > 0)

def ratio(numerator, reference):
    numerator, reference = np.broadcast_arrays(np.asarray(numerator, dtype=float), np.asarray(reference, dtype=float))
    result = np.full(numerator.shape, np.nan)
    valid = ~np.isnan(numerator) & ~np.isnan(reference)
    np.divide(numerator, reference, out=result, where=valid & (reference != 0))
    result[valid & (reference == 0) & (numerator > 0)] = np.inf
    return result

def interval(values):
    values = np.asarray(values, dtype=float)
    bad = int(np.isnan(values).sum())
    inf = int(np.isposinf(values).sum())
    if bad:
        return (np.nan, np.nan, bad, inf, 'unavailable_undefined_draw')
    ordered = np.sort(values)

    def percentile(q):
        position = (len(ordered) - 1) * q
        lo = int(np.floor(position))
        hi = int(np.ceil(position))
        w = position - lo
        if lo == hi or ordered[lo] == ordered[hi]:
            return ordered[lo]
        if np.isposinf(ordered[hi]):
            return np.inf
        return ordered[lo] * (1 - w) + ordered[hi] * w
    return (percentile(0.025), percentile(0.975), 0, inf, 'extended_percentile' if inf else 'finite_percentile')

def endpoint_summary(members, design, available):
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
    scopes = [('run', sample, [sample]) for sample in names] + [('organism', sp, [x for x in names if present[x] == sp]) for sp in BITS_SPECIES] + [('pooled', 'all', names)]
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
BITS_SPECIES = ['H', 'Y', 'E']

def bootstrap(summary, a, p, pairs, keys, design, out, draws=2000):
    """One shared recipient/key draw per replicate across all filters and endpoints."""
    if draws != 2000:
        raise ValueError('Frozen draw count is 2000')
    singles = design[design.kind.eq('single')]
    names = singles['sample'].tolist()
    strata = [np.flatnonzero(singles.present_species.to_numpy() == sp) for sp in BITS_SPECIES]
    ri = pairs.run_index.to_numpy(dtype=int)
    ki = pairs.key_index.to_numpy(dtype=int)
    rng = np.random.Generator(np.random.PCG64(SEED))
    nkeys = len(keys)
    ledger = {scheme: np.zeros((2, len(summary), draws), dtype=np.int64) for scheme in ['recipient', 'precursor', 'two_way']}
    for start in range(0, draws, 16):
        width = min(16, draws - start)
        rw = np.zeros((18, width), dtype=np.int64)
        kw = np.zeros((nkeys, width), dtype=np.int64)
        for column in range(width):
            for group in strata:
                rw[:, column] += np.bincount(rng.choice(group, len(group), replace=True), minlength=18)
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

def contrasts(summary, ledger, out):
    rows = []
    drawrows = []
    rawdraws = []
    dims = ['endpoint', 'mask', 'scope', 'group', 'tercile']
    for _, group in summary.groupby(dims, sort=False):
        if group.endpoint.iloc[0] == 'T1':
            indexed = group.set_index('class')
            pairs = [('T', 'R'), ('L', 'R')]
        else:
            indexed = group.set_index('support_bin')
            pairs = [(BINS[-1], BINS[1])]
        for numerator, reference in pairs:
            if numerator not in indexed.index or reference not in indexed.index:
                continue
            n = indexed.loc[numerator]
            r = indexed.loc[reference]
            base = {x: n[x] for x in dims}
            pointn = float(n.absent_fraction)
            pointr = float(r.absent_fraction)
            base.update({'numerator': numerator, 'reference': reference, 'numerator_A': n.A, 'numerator_n': n.n, 'reference_A': r.A, 'reference_n': r.n, 'fraction_difference': pointn - pointr, 'fraction_ratio': float(ratio(pointn, pointr))})
            rows.append(base)
            for scheme, counts in ledger.items():
                nf = fraction(*counts[:, int(n.endpoint_id), :])
                rf = fraction(*counts[:, int(r.endpoint_id), :])
                for statistic, values in [('difference', nf - rf), ('ratio', ratio(nf, rf))]:
                    low, high, bad, inf, state = interval(values)
                    drawrows.append({**base, 'contrast_draw_id': len(rawdraws), 'scheme': scheme, 'statistic': statistic, 'low95': low, 'high95': high, 'undefined_draws': bad, 'infinite_draws': inf, 'interval_status': state, 'draws': 2000})
                    rawdraws.append(values)
    table(out / 'contrasts.tsv', pd.DataFrame(rows))
    table(out / 'contrast_intervals.tsv', pd.DataFrame(drawrows))
    np.savez_compressed(out / 'contrast_draws.npz', values=np.asarray(rawdraws))
