"""Supplementary Data 1 (tables D1 to D11) from the pipeline outputs, plus the v5.1 tables behind the short
Supplement tables S3, S5, S6 and S10.

Each function returns one sheet as a list of blocks (title, DataFrame). Values are kept at full precision; the
printed rounding of the published workbook is applied only when comparing (scripts/supplementary_data.py).
Column names and row labels follow the published workbook. Labels of the v5.1 vocabulary are kept as published
(second-pass-only = added, first-pass-only = lost).
"""
from pathlib import Path

import numpy as np
import pandas as pd

CONF = ['on', 'off', 'standard']
LONG = {'on': 'original', 'off': 'no_refine', 'standard': 'standard_library'}
CLASS_LABEL = {'first_pass': 'First pass', 'final': 'final', 'retained': 'Retained', 'second_pass_only': 'Second-pass-only',
               'lost': 'First-pass-only', 'lost_library_present': 'First-pass-only, in library',
               'lost_library_excluded': 'First-pass-only, outside library'}
LEVEL_LABEL = {'precursor': 'Precursor', 'peptide': 'Peptide', 'protein_group': 'Protein group'}
SPECIES = {'H': 'Human', 'Y': 'Yeast', 'E': 'E. coli'}


def _read(path, **kw):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kw)


def _pair(a, b):
    return f'{int(a)}/{int(b)}'


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


def _num(x):
    x = _f(x)
    return 'unavailable' if not np.isfinite(x) else repr(x)


def _to(lo, hi):
    return f'{_num(lo)} to {_num(hi)}'


def _br(point, lo, hi):
    return f'{_num(point)} [{_num(lo)}, {_num(hi)}]'


class Sources:
    def __init__(self, root):
        self.root = Path(root)
        self.sc = self.root / 'single-cell'
        self.pg = self.root / 'pg250'

    def acc(self, conf, name):
        return _read(self.sc / conf / 'accounting' / name)

    def iv(self, conf, name):
        return _read(self.sc / conf / 'intervals' / name)






def _census(s):
    parts = [_read(s.sc / 'empirical-census-on.tsv', dtype={'charge': str})]
    predicted = s.sc / 'predicted-census.tsv'
    if predicted.exists():
        parts.append(_read(predicted, dtype={'charge': str}))
    return pd.concat(parts, ignore_index=True)


def s9a(s):
    c = _census(s)
    c = c[(c.unit == 'precursor') & (c.classification == 'I_as_L_exclusive')]
    rows = [{'Library': r.library, 'Charge': r.charge, 'Total keys': r.total_candidates, 'Human': r.human_only, 'Yeast': r.yeast_only,
             'Shared': r.shared_accessions, 'I/L shared': r.shared_IL, 'r': r.r} for r in c.itertuples()]
    return [('', pd.DataFrame(rows))]












def _endpoint(s, folder):
    counts = _read(folder / 'endpoint_counts.tsv', keep_default_na=False)
    iv = _read(folder / 'fraction_intervals.tsv', keep_default_na=False)
    two = iv[iv.scheme == 'two_way'].set_index('endpoint_id')
    return counts, two










def _scales(s, repeat):
    return _read(s.sc / ('comparison-' + repeat) / 'class_species_scales_comparison.tsv')






def s49(s):
    t = _read(s.sc / 'comparison-off' / 'empirical_library_opportunities.tsv')
    names = {'original_default': 'Refinement on', 'disabled_refinement': 'Refinement off'}
    rows = [{'Configuration': names[r.configuration], 'Candidate keys': r.total_candidates, 'Human only': r.human_only,
             'Yeast only': r.yeast_only, 'Shared': r.shared, 'I/L-shared': r.shared_IL, 'Unknown': r.unknown, 'Yeast/human ratio': r.r}
            for _, r in t.iterrows()]
    return [('', pd.DataFrame(rows))]




def s72(s):
    sources = {'on': (_scales(s, 'off'), 'original_'), 'off': (_scales(s, 'off'), 'norefine_'), 'standard': (_scales(s, 'standard'), 'norefine_')}
    rows = []
    for cls in ['retained', 'second_pass_only', 'lost', 'lost_library_present', 'lost_library_excluded']:
        for conf in CONF:
            table, prefix = sources[conf]
            r = table[(table['filter'] == 'F3') & (table['class'] == cls)].iloc[0]
            iv = s.iv(conf, 'class_intervals.tsv')
            i = iv[(iv['class'] == cls) & (iv.scheme == 'two_way')].iloc[0]
            rows.append({'Observation class / setting': CLASS_LABEL[cls] + ' / ' + conf,
                         'Yeast/n': _pair(r[prefix + 'yeast'], r[prefix + 'species_exclusive_denominator']),
                         'Yeast,%': 100 * r[prefix + 'observed_yeast_fraction'], 'Two-way 95% interval,%': _to(100 * i.low95, 100 * i.high95),
                         'Ratio / retained': r[prefix + 'r_free_ratio_vs_retained'],
                         'Historical scale,%': 100 * r[prefix + 'conditional_historical'],
                         'Empirical scale,%': 100 * r[prefix + 'conditional_configuration_empirical']})
    return [('', pd.DataFrame(rows))]


def _long(s, name, keys, folder='accounting', metrics=None):
    frames = {}
    for conf in CONF:
        t = (s.acc if folder == 'accounting' else s.iv)(conf, name)
        value_columns = [c for c in t.columns if c not in keys and c != 'primary'] if metrics is None else metrics
        frames[LONG[conf]] = t.melt(id_vars=keys, value_vars=value_columns, var_name='metric', value_name='value')
    out = frames['original'].rename(columns={'value': 'original'})
    for name2 in ['no_refine', 'standard_library']:
        out[name2] = frames[name2]['value'].to_numpy()
    out['standard_library_status'] = 'complete'
    return [('', out)]


def s61(s):
    return _long(s, 'accounting.tsv', ['filter'])


def s62(s):
    return _long(s, 'class_intervals.tsv', ['class', 'scheme'], 'intervals',
                 ['point', 'low95', 'high95', 'draws', 'undefined_draws', 'infinite_draws', 'interval_status', 'yeast', 'species_exclusive_denominator', 'precursor_blocks'])


def s63(s):
    return _long(s, 'class_ratio_intervals.tsv', ['contrast', 'scheme'], 'intervals',
                 ['point', 'low95', 'high95', 'draws', 'undefined_draws', 'infinite_draws', 'interval_status', 'numerator_yeast', 'numerator_denominator', 'reference_yeast', 'reference_denominator', 'precursor_blocks'])


def s65(s):
    return _long(s, 'represented_coverage.tsv', ['filter', 'level'])




def s71(s):
    ent = {c: s.acc(c, 'entity_species.tsv') for c in CONF}
    cols = {'on': 'Refinement on', 'off': 'Refinement off', 'standard': 'Standard library, refinement off'}
    rows = []
    for level in ['precursor', 'peptide', 'protein_group']:
        get = {c: {st: ent[c][(ent[c]['filter'] == 'F3') & (ent[c].level == level) & (ent[c].state == st)].iloc[0]
                   for st in ['wholly_absent', 'still_represented']} for c in CONF}
        for label, fn in [('Wholly absent', lambda g: int(g['wholly_absent'].observations)),
                          ('Human-only absent', lambda g: int(g['wholly_absent'].human)),
                          ('Absent yeast/n', lambda g: _pair(g['wholly_absent'].yeast, g['wholly_absent'].species_exclusive_denominator)),
                          ('Absent yeast,%', lambda g: 100 * g['wholly_absent'].observed_yeast_fraction),
                          ('Still represented yeast/n', lambda g: _pair(g['still_represented'].yeast, g['still_represented'].species_exclusive_denominator)),
                          ('Still yeast,%', lambda g: 100 * g['still_represented'].observed_yeast_fraction),
                          ('Absent/still fraction ratio', lambda g: g['wholly_absent'].absent_vs_still_yeast_fraction_ratio)]:
            rows.append({'Entity': LEVEL_LABEL[level], 'Measure': label, **{cols[c]: fn(get[c]) for c in CONF}})
    return [('', pd.DataFrame(rows))]




def _pool_alias(name):
    if name.startswith('homogeneous_6_'):
        return 'H6-' + name[-1]
    if name.startswith('mixed_6_'):
        return 'M6-' + name[-1]
    if name.startswith('homogeneous_16_'):
        return 'H16-' + name.split('_')[-1]
    if name.startswith('mixed_16_'):
        return 'M16-' + name[-1]
    return name


def s1b(s):
    m = _read(s.sc / 'pools' / 'pool_membership.tsv')
    rules = {('H', 6): 'one cell line, first six planned ranks', ('M', 6): 'four lines with counts 2+2+1+1; fixed planned membership',
             ('H', 16): 'one cell line, planned ranks 1 to 16', ('M', 16): 'four cells per line; zero-based within-line rank modulo 4 equals pool index minus 1'}
    rows = []
    for pool, g in m.groupby('pool', sort=False):
        arm = 'H' if pool.startswith('homogeneous') else 'M'
        size = len(g)
        counts = g.line.value_counts().sort_index()
        rows.append({'Pool': _pool_alias(pool), 'Size': size, 'Arm': arm, 'Composition': '; '.join(f'{k}:{v}' for k, v in counts.items()),
                     'Rule': rules[arm, size]})
    return [('', pd.DataFrame(rows))]


def s2(s):
    t = _read(s.sc / 'pools' / 'A_composition_uncertainty.tsv')
    names = {'conditional_count_model': 'Original count model', 'recipient': 'Recipient-paired', 'precursor': 'Precursor blocks', 'two_way': 'Two-way'}
    rows = [{'Size': int(r.cohort.split('_')[-1]), 'Method': names[r.scheme], 'Ratio': r.point, '95% interval': _to(r.low95, r.high95),
             'Recipients': int(r.recipients), 'Precursor blocks': int(r.precursor_blocks)} for _, r in t.iterrows()]
    return [('', pd.DataFrame(rows))]






def s4(s):
    p = _read(s.sc / 'pools' / 'A_per_pool.tsv')
    rows = [{'Size': int(r.pool_size), 'Arm': r.composition, 'Pool': _pool_alias(r.pool), 'Recipients': int(r.recipients),
             'Second-pass-only observations': int(r.events), 'Yeast': int(r.yeast), 'Yeast/second-pass-only observation,%': 100 * r.observed_yeast_per_transfer}
            for _, r in p.iterrows()]
    return [('', pd.DataFrame(rows))]


T1_CLASS = {'D': 'First pass', 'F': 'Final', 'R': 'Retained', 'T': 'Second-pass-only', 'L': 'First-pass-only',
            'L_in_library': 'First-pass-only, in library', 'L_outside_library': 'First-pass-only, outside library'}


def s51(s, design):
    counts, two = _endpoint(s, s.pg / 'primary' / 'export')
    order = {name: i + 1 for i, name in enumerate(design['sample'])}
    sel = counts[(counts.endpoint == 'T1') & (counts['mask'] == 'F3') & (counts.scope == 'run')]
    rows = []
    for _, r in sel.iterrows():
        i = two.loc[int(r.endpoint_id)]
        rows.append({'Run': f'Run {order[r["group"]]}', 'Class': T1_CLASS[r['class']], 'A': int(r.A), 'P': int(r.P), 'n': int(r.n),
                     'Shared': int(r.excluded_shared), 'I/L-shared': int(r.excluded_shared_IL), 'Unknown': int(r.excluded_unknown),
                     'Quarantined': int(r.excluded_quarantined),
                     'Observed wrong-organism fraction [two-way 95% interval]': _br(r.absent_fraction, i.low95, i.high95)})
    return [('', pd.DataFrame(rows))]


def s73(s):
    fold = s.pg / 'primary' / 'fold-change'
    t = _read(fold / 't3_summary.tsv', keep_default_na=False)
    c = _read(fold / 't3_coverage.tsv', keep_default_na=False)
    t = t[(t['mask'] == 'F3') & (t.scheme == 'stage_specific') & (t.stage == 'final') & (t.transition == 'F')]
    rows = []
    for org in ['H', 'Y', 'E', 'ALL']:
        for threshold, tlabel in [(0, 'Unfiltered'), (3, 'At least 3/6'), (5, 'At least 5/6')]:
            for unit in ['rows', 'events']:
                if org == 'ALL':
                    scope = 'pooled_organism_key' if unit == 'rows' else 'pooled_events'
                    r = t[(t.organism == org) & (t.threshold == threshold) & (t.unit == unit) & (t.scope == scope)].iloc[0]
                    label = f'{tlabel} organism-key rows' if unit == 'rows' else f'{tlabel} events'
                    kept = 'not supplied as pooled coverage'
                else:
                    r = t[(t.organism == org) & (t.threshold == threshold) & (t.unit == unit) & (t.scope == 'group')].iloc[0]
                    label = f'{tlabel} {unit}'
                    cv = c[(c['mask'] == 'F3') & (c.organism == org) & (c.scheme == 'stage_specific') & (c.stage == 'final') & (c.threshold == threshold) & (c.transition == 'F')].iloc[0]
                    kept = _pair(cv.present_rows_retained, cv.present_rows_baseline) if unit == 'rows' else _pair(cv.present_events_retained, cv.present_events_baseline)
                rows.append({'Intended organism': 'Pooled' if org == 'ALL' else SPECIES[org], 'Population / unit': label,
                             'A': int(r.A), 'P': int(r.P), 'n': int(r.n),
                             'Observed wrong-organism fraction [two-way 95% interval]': _br(r.false_fraction, r.two_way_ci_low, r.two_way_ci_high),
                             'Intended-species retained / baseline': kept,
                             'Precursor-only 95% interval (rows; post hoc)': _to(r.precursor_ci_low, r.precursor_ci_high) if unit == 'rows' else 'not applicable'})
    return [('', pd.DataFrame(rows))]


def s60(s):
    a, b = [], []
    for org in ['H', 'Y', 'E']:
        folder = s.pg / ('within-' + org)
        points = _read(folder / 'contrasts.tsv')
        limits = _read(folder / 'contrast_intervals.tsv')
        counts = _read(folder / 'endpoint_counts.tsv', keep_default_na=False)
        for mask in ['F3', 'F0', 'F4']:
            for endpoint, num, ref, label, tercile in [('T1', 'T', 'R', 'T/R', 'all'), ('T1', 'L', 'R', 'L/R', 'all'),
                                                         ('T2', 'positive_60_100', 'positive_0_20', 'high/low support, high quantity', 'high')]:
                pc = (points.endpoint == endpoint) & (points['mask'] == mask) & (points.scope == 'organism') & (points.numerator == num) & (points.reference == ref) & (points.tercile.astype(str) == tercile)
                p = points[pc].iloc[0]
                lc = (limits.endpoint == endpoint) & (limits['mask'] == mask) & (limits.scope == 'organism') & (limits.numerator == num) & (limits.reference == ref) & (limits.tercile.astype(str) == tercile) & (limits.scheme == 'two_way') & (limits.statistic == 'ratio')
                lim = limits[lc].iloc[0]
                a.append({'Cohort': org, 'Mask': mask, 'Contrast': label,
                          'Numerator A/P': _pair(p.numerator_A, p.numerator_n - p.numerator_A), 'Reference A/P': _pair(p.reference_A, p.reference_n - p.reference_A),
                          'Fraction ratio': p.fraction_ratio, 'Fraction difference': p.fraction_difference,
                          '95% two-way interval': _to(lim.low95, lim.high95) if lim.interval_status == 'finite_percentile' else 'unavailable',
                          'Undefined draws': int(lim.undefined_draws), 'Infinite draws': int(lim.infinite_draws)})
            for cls in ['D', 'F', 'R', 'T', 'L']:
                r = counts[(counts.endpoint == 'T1') & (counts['mask'] == mask) & (counts.scope == 'organism') & (counts['class'] == cls)].iloc[0]
                b.append({'Cohort': org, 'Mask': mask, 'Class': cls, 'events': int(r.events), 'A': int(r.A), 'P': int(r.P), 'n': int(r.n),
                          'absent fraction': float(r.absent_fraction)})
    return [('Block A', pd.DataFrame(a)), ('Block B', pd.DataFrame(b))]


def s56(s):
    q = _read(s.pg / 'native-qa' / 'native-run-qa.tsv', keep_default_na=False)
    q = q[q.kind == 'blank']
    rows = [{'Pass': r['pass'], 'Ordinal': int(r.manifest_ordinal), 'Native sample': r['sample'], 'Identified precursors (QA only)': int(r['Precursors.Identified']),
             'MS1': r.search_ms1_ppm, 'MS2': r.search_ms2_ppm, 'Calibration failures': int(r.calibration_failures)} for _, r in q.iterrows()]
    return [('', pd.DataFrame(rows))]


CONTEXT = {'D': 'First pass', 'F': 'Final', 'R': 'Retained', 'T': 'Second-pass-only'}


def _fc(s, study):
    f = _read(s.pg / study / 'fold-change' / 'fc_summary.tsv', keep_default_na=False)
    return f[(f['mask'] == 'F3') & (f.quantity == 'Precursor.Normalised')]


def s53(s):
    f = _fc(s, 'primary')
    rows = []
    for sp in ['Y', 'E']:
        for ctx in 'DFRT':
            r = f[(f.species == sp) & (f.context == ctx) & (f.population == 'own_context')].iloc[0]
            rows.append({'Organism': SPECIES[sp], 'Context': CONTEXT[ctx], 'Eligible keys': int(r.eligible_keys), 'Universe keys': int(r.universe_keys),
                         'Expected': float(r.expected_log2fc), 'Observed median': float(r.median_observed_log2fc), 'Median bias': float(r.median_bias),
                         'MAE': float(r.MAE), 'RMSE': float(r.RMSE)})
    return [('', pd.DataFrame(rows))]


def s54(s):
    f = _fc(s, 'primary')
    rows = []
    for sp in ['E', 'Y']:
        for pop, ctxs in [('D_F_intersection', 'DF'), ('R_T_intersection', 'RT')]:
            for ctx in ctxs:
                r = f[(f.species == sp) & (f.context == ctx) & (f.population == pop)].iloc[0]
                rows.append({'Organism': SPECIES[sp], 'Context': CONTEXT[ctx], 'Paired population': pop, 'Eligible keys': int(r.eligible_keys),
                             'Universe keys': int(r.universe_keys), 'MAE': float(r.MAE), 'RMSE': float(r.RMSE)})
    return [('', pd.DataFrame(rows))]


def s58(s):
    a, b = _fc(s, 'primary'), _fc(s, 'replication')
    rows = []
    for sp in ['Y', 'E']:
        for ctx in 'DFRT':
            x = a[(a.species == sp) & (a.context == ctx) & (a.population == 'own_context')].iloc[0]
            y = b[(b.species == sp) & (b.context == ctx) & (b.population == 'own_context')].iloc[0]
            rows.append({'Species': SPECIES[sp], 'Context': CONTEXT[ctx], 'Jan 16 median': float(x.median_observed_log2fc), 'Jan 16 bias': float(x.median_bias),
                         'Jan 16 MAE': float(x.MAE), 'Jan 16 keys': int(x.eligible_keys), 'Jan 9 median': float(y.median_observed_log2fc),
                         'Jan 9 bias': float(y.median_bias), 'Jan 9 MAE': float(y.MAE), 'Jan 9 keys': int(y.eligible_keys)})
    block = []
    for sp in ['E', 'Y']:
        for pop, ctxs in [('D_F_intersection', 'DF'), ('R_T_intersection', 'RT')]:
            for ctx in ctxs:
                y = b[(b.species == sp) & (b.context == ctx) & (b.population == pop)].iloc[0]
                block.append({'Species': SPECIES[sp], 'Paired population': pop, 'Context': CONTEXT[ctx], 'Eligible keys': int(y.eligible_keys),
                              'Median log2 fold change': float(y.median_observed_log2fc), 'Median bias': float(y.median_bias), 'MAE': float(y.MAE)})
    return [('', pd.DataFrame(rows)), ('Continued', pd.DataFrame(block))]


METRIC = {'median_delta': 'Median delta', 'fraction_toward_zero': 'Toward zero', 'MAE_change': 'MAE change', 'MAE_R': 'MAE retained',
          'MAE_F': 'MAE final', 'fraction_tie': 'Tie', 'fraction_away_from_zero': 'Away from zero', 'missed_twofold': 'Missed two-fold',
          'wrong_sign': 'Wrong sign', 'within_half': 'Within 0.5', 'median_bias_R': 'Median bias R', 'median_bias_T': 'Median bias T',
          'median_bias_T_minus_R': 'Bias T minus R'}
FRACTIONS = {'fraction_toward_zero', 'fraction_tie', 'fraction_away_from_zero', 'missed_twofold', 'wrong_sign', 'within_half'}


def s59(s):
    t = _read(s.pg / 'posthoc-fc' / 'summary.tsv', keep_default_na=False)
    batch = {'jan16': '2025-01-16', 'jan09': '2025-01-09'}
    rows = []
    for _, r in t.iterrows():
        pct = r.metric in FRACTIONS
        scale = 100 if pct else 1
        unit = '%' if pct else ''
        point = float(r.point) * scale
        lo, hi = float(r.ci_low) * scale, float(r.ci_high) * scale
        pop = r.population.split('|')[-1]
        label = {'P1': 'P1: paired R/F'}.get(pop, pop.replace('_', ': ', 1))
        dash = lambda v: '-' if v in ('', 'NA') else int(float(v))
        rows.append({'Batch': batch[r.batch], 'Mask': r['mask'], 'Species': SPECIES[r.species],
                     'Quantity': 'Normalised' if r.quantity == 'Precursor.Normalised' else 'Raw', 'Endpoint and population': label,
                     'Metric': METRIC[r.metric], 'Point': f'{point!r}{unit}' if pct else point, '95% interval': f'[{lo!r}{unit}, {hi!r}{unit}]',
                     'Undefined draws': int(r.undefined_draws), 'Status': r.status, 'Universe': int(r.universe_keys),
                     'Eligible': dash(r.get('eligible_keys', '')), 'Numerator': dash(r.get('numerator', '')), 'Union': dash(r.get('union_keys', '')),
                     'R': dash(r.get('R_keys', '')), 'T': dash(r.get('T_keys', '')), 'Shared': dash(r.get('shared_keys', ''))})
    return [('', pd.DataFrame(rows))]


def s64(s):
    """Class composition and scales for every filter, class and configuration (v5.1 Table S64)."""
    off, std = _scales(s, 'off'), _scales(s, 'standard').set_index(['filter', 'class'])
    metrics = ['observations', 'human', 'yeast', 'species_exclusive_denominator', 'observed_yeast_fraction', 'r_free_ratio_vs_retained',
               'historical_r', 'configuration_empirical_r', 'conditional_historical', 'conditional_configuration_empirical',
               'same_cell_group_in_final_fraction', 'same_cell_group_other_peptide_in_final_fraction']
    rows = []
    for _, r in off.iterrows():
        t = std.loc[(r['filter'], r['class'])]
        for m in metrics:
            rows.append({'filter': r['filter'], 'class': r['class'], 'metric': m, 'original': r['original_' + m],
                         'no_refine': r['norefine_' + m], 'standard_library': t['norefine_' + m], 'standard_library_status': 'complete'})
    return [('', pd.DataFrame(rows))]


# Sheets of the published workbook (D1 to D11), each the v5.1 tables it holds, in block order, and the number of
# leading label columns that identify a row in each block.
SHEETS = {'D1': [('s61', 2)], 'D2': [('s64', 3)], 'D3': [('s62', 3), ('s63', 3)], 'D4': [('s65', 3)],
          'D5': [('s9a', 2), ('s49', 1)], 'D6': [('s1b', 1), ('s4', 3)], 'D7': [('s51', 2)], 'D8': [('s73', 2)],
          'D9': [('s60', 3)], 'D10': [('s58', [2, 3]), ('s53', 2), ('s54', 3)], 'D11': [('s59', 6)]}
# v5.1 tables cited by the short Supplement but not held in the workbook; written as extra sheets.
EXTRA = {'Table S3 (v5.1 S72)': [('s72', 1)], 'Table S5 (v5.1 S71)': [('s71', 2)], 'Table S6 (v5.1 S2)': [('s2', 2)],
         'Table S10 (v5.1 S56)': [('s56', 2)]}
KEYS = {}


def build(root, design_path):
    s = Sources(root)
    design = pd.read_csv(design_path, sep='\t', keep_default_na=False)
    sheets = {}
    for name, parts in {**SHEETS, **EXTRA}.items():
        blocks, keys = [], []
        try:
            for function, nkeys in parts:
                fn = globals()[function]
                produced = fn(s, design) if function == 's51' else fn(s)
                blocks += produced
                keys += nkeys if isinstance(nkeys, list) else [nkeys] * len(produced)
            sheets[name] = blocks
            KEYS[name] = keys
        except Exception as error:  # a missing step output leaves its sheet unproduced, reported as such
            sheets[name] = error
    return sheets
