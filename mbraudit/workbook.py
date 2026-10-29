"""Supplementary Data 1 (tables D1 to D12) from the pipeline outputs, plus the v5.1 tables behind the short
Supplement tables S3, S5, S6 and S10. D12 holds the source data of main-text Figure 2 (Tables 3 and 4 of the earlier
manuscript, m3 and m4).

Each function returns one sheet as a list of blocks (title, DataFrame). Values are kept at full precision (the
long tables D2 to D5 as the exact text of the output tables). publish() writes the published workbook: a README
sheet and one sheet per table with its title, description and block labels (PUBLICATION), each value written with
the digits of the published tables (FORMATS). Column names and row labels follow the published workbook. Labels of
the v5.1 vocabulary are kept as published (second-pass-only = added, first-pass-only = lost).
"""
import math
import re
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


def _text(path):
    """An output table as its exact text (no number parsing)."""
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


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

    def acc(self, conf, name, text=False):
        return (_text if text else _read)(self.sc / conf / 'accounting' / name)

    def iv(self, conf, name, text=False):
        return (_text if text else _read)(self.sc / conf / 'intervals' / name)






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










def _scales(s, repeat, text=False):
    return (_text if text else _read)(s.sc / ('comparison-' + repeat) / 'class_species_scales_comparison.tsv')






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
            r = table[(table['filter'] == 'F1') & (table['class'] == cls)].iloc[0]
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
        t = (s.acc if folder == 'accounting' else s.iv)(conf, name, text=True)
        value_columns = [c for c in t.columns if c not in keys and c != 'primary'] if metrics is None else metrics
        long = t.melt(id_vars=keys, value_vars=value_columns, var_name='metric', value_name='value', ignore_index=False)
        # one source row after the other, each with its metrics in the listed order (the published layout)
        long['_row'], long['_metric'] = long.index, long.metric.map({m: i for i, m in enumerate(value_columns)})
        frames[LONG[conf]] = long.sort_values(['_row', '_metric'], kind='stable').drop(columns=['_row', '_metric']).reset_index(drop=True)
    out = frames['original'].rename(columns={'value': 'original'})
    for name2 in ['no_refine', 'standard_library']:
        out[name2] = frames[name2]['value'].to_numpy()
    out['standard_library_status'] = 'complete'
    return [('', out)]


def s61(s):
    return _long(s, 'accounting.tsv', ['filter'], metrics=[
        'first_observations', 'final_observations', 'retained_observations', 'second_pass_only_observations', 'lost_observations',
        'lost_library_present', 'lost_library_excluded', 'net_observation_change_fraction', 'second_pass_only_with_F3_firstpass',
        'lost_with_F3_final'])


def s62(s):
    return _long(s, 'class_intervals.tsv', ['class', 'scheme'], 'intervals',
                 ['point', 'low95', 'high95', 'draws', 'undefined_draws', 'infinite_draws', 'interval_status', 'yeast', 'species_exclusive_denominator', 'precursor_blocks'])


def s63(s):
    return _long(s, 'class_ratio_intervals.tsv', ['contrast', 'scheme'], 'intervals',
                 ['point', 'low95', 'high95', 'draws', 'undefined_draws', 'infinite_draws', 'interval_status', 'numerator_yeast', 'numerator_denominator', 'reference_yeast', 'reference_denominator', 'precursor_blocks'])


def s65(s):
    return _long(s, 'represented_coverage.tsv', ['filter', 'level'], metrics=[
        'first_represented_entities', 'final_represented_entities', 'common_entities', 'first_all_row_completeness',
        'final_all_row_completeness', 'first_common_row_completeness', 'final_common_row_completeness',
        'first_positive_common_entity_cells', 'final_positive_common_entity_cells', 'common_row_denominator'])




def s71(s):
    ent = {c: s.acc(c, 'entity_species.tsv') for c in CONF}
    cols = {'on': 'Refinement on', 'off': 'Refinement off', 'standard': 'Standard library, refinement off'}
    rows = []
    for level in ['precursor', 'peptide', 'protein_group']:
        get = {c: {st: ent[c][(ent[c]['filter'] == 'F1') & (ent[c].level == level) & (ent[c].state == st)].iloc[0]
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
    frame = pd.DataFrame(rows).sort_values(['Size', 'Arm', 'Pool'], kind='stable').reset_index(drop=True)
    return [('', frame)]


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
    sel = counts[(counts.endpoint == 'T1') & (counts['mask'] == 'F1') & (counts.scope == 'run')]
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
    t = t[(t['mask'] == 'F1') & (t.scheme == 'stage_specific') & (t.stage == 'final') & (t.transition == 'F')]
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
                    cv = c[(c['mask'] == 'F1') & (c.organism == org) & (c.scheme == 'stage_specific') & (c.stage == 'final') & (c.threshold == threshold) & (c.transition == 'F')].iloc[0]
                    kept = _pair(cv.present_rows_retained, cv.present_rows_baseline) if unit == 'rows' else _pair(cv.present_events_retained, cv.present_events_baseline)
                rows.append({'Intended organism': 'Pooled' if org == 'ALL' else SPECIES[org], 'Population / unit': label,
                             'A': int(r.A), 'P': int(r.P), 'n': int(r.n),
                             'Observed wrong-organism fraction [two-way 95% interval]': _br(r.false_fraction, r.two_way_ci_low, r.two_way_ci_high),
                             'Intended-species retained / baseline': kept,
                             'Precursor-only 95% interval (rows; post hoc)': _to(r.precursor_ci_low, r.precursor_ci_high) if unit == 'rows' else 'not applicable'})
    return [('', pd.DataFrame(rows))]


def reading(row):
    """The frozen reading rule of a within-organism contrast (carried from the report of that analysis)."""
    try:
        ratio = float(row['fraction_ratio'])
    except (ValueError, TypeError):
        ratio = math.nan
    if not math.isfinite(ratio):
        return 'Point ratio unavailable or infinite; no finite direction claim.'
    if row['endpoint'] == 'T2':
        point = 'Below one: agrees directionally with low-support enrichment.' if ratio < 1 else ('Above one: reverses the low-support direction.' if ratio > 1 else 'Point ratio equals one.')
    elif row['numerator'] == 'T':
        point = 'Higher observed wrong-organism fraction in second-pass-only events.' if ratio > 1 else ('Lower observed wrong-organism fraction in second-pass-only events.' if ratio < 1 else 'Point fractions equal.')
    else:
        point = 'Above one: aligns directionally with Astral L/R.' if ratio > 1 else ('Below one: aligns with the earlier mixed-species reversal.' if ratio < 1 else 'Point fractions equal.')
    try:
        lo, hi = float(row['low95']), float(row['high95'])
    except (ValueError, TypeError):
        lo = hi = math.nan
    if math.isfinite(lo) and math.isfinite(hi):
        point += ' Interval includes no difference.' if lo <= 1 <= hi else ' Finite interval excludes one.'
    else:
        point += ' Interval unavailable or infinite; finite interval interpretation unavailable.'
    return point


def s60(s):
    a, b = [], []
    for org in ['H', 'Y', 'E']:
        folder = s.pg / ('within-' + org)
        points = _read(folder / 'contrasts.tsv')
        limits = _read(folder / 'contrast_intervals.tsv')
        counts = _read(folder / 'endpoint_counts.tsv', keep_default_na=False)
        for mask in ['F1', 'F2', 'F3']:
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
                          'Undefined draws': int(lim.undefined_draws), 'Infinite draws': int(lim.infinite_draws),
                          'Frozen-rule reading': reading(lim), 'Source line': f'contrast_intervals.tsv:{lim.name + 2}'})
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
    return f[(f['mask'] == 'F1') & (f.quantity == 'Precursor.Normalised')]


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
    for pop, ctxs in [('D_F_intersection', 'DF'), ('R_T_intersection', 'RT')]:
        for sp in ['E', 'Y']:
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
                     'R': dash(r.get('R_keys', '')), 'T': dash(r.get('T_keys', '')), 'Shared': dash(r.get('shared_keys', '')),
                     'Source line': r.name + 2})
    # Blocks of the published table: the primary rows of P1, P2 and P3 in the order of summary.tsv, then every
    # sensitivity row by endpoint, batch, filter (F1, F2, F3), quantity and species.
    from mbraudit.filters import FILTERS as MASKS
    from mbraudit.pg250.posthoc_fc import QUANTITIES
    frame = pd.DataFrame(rows)
    endpoint = t.population.str.split('|').str[-1].str[:2]
    primary = (t['display'] == 'primary').to_numpy()
    order = pd.DataFrame({'endpoint': endpoint, 'batch': t.batch.map({'jan16': 0, 'jan09': 1}),
                          'mask': t['mask'].map({m: i for i, m in enumerate(MASKS)}),
                          'quantity': t.quantity.map({q: i for i, q in enumerate(QUANTITIES)}), 'species': t.species.map({'Y': 0, 'E': 1})})
    sensitivity = order[~primary].sort_values(list(order.columns), kind='stable').index
    endpoint = endpoint.to_numpy()
    return [('', frame[primary & (endpoint == 'P1')].reset_index(drop=True)),
            ('P2 primary', frame[primary & (endpoint == 'P2')].reset_index(drop=True)),
            ('P3 primary', frame[primary & (endpoint == 'P3')].reset_index(drop=True)),
            ('Sensitivity', frame.loc[sensitivity].reset_index(drop=True))]


def s64(s):
    """Class composition and scales for every filter, class and configuration (v5.1 Table S64)."""
    off, std = _scales(s, 'off', text=True), _scales(s, 'standard', text=True).set_index(['filter', 'class'])
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


def m3(s):
    """Table 3 of the earlier manuscript, now the source data of Figure 2A (D12 Block A): all 54 runs, F1."""
    from mbraudit.tables import wrong_organism
    rows = [{'Intended organism': r['organism'], 'Class': r['class'], 'A': r['A'], 'P': r['P'], 'n': r['n'],
             'Observed wrong-organism fraction [two-way 95% interval]': _br(r['fraction'], r['low95'], r['high95'])}
            for r in wrong_organism(s.pg / 'primary' / 'export')]
    return [('', pd.DataFrame(rows))]


def m4(s):
    """Table 4 of the earlier manuscript, now the source data of Figure 2B (D12 Block B): the human-only runs
    analysed with all 54 runs and alone, F1."""
    from mbraudit.tables import human_only_passes
    passes, _ = human_only_passes(s.pg / 'primary' / 'export', s.pg / 'within-H')
    rows = [{'Class': r['class'], 'With all 54 runs: A/n': f"{r['with_A']:,}/{r['with_n']:,}", 'With all 54 runs: %': r['with_percent'],
             'Human-only runs alone: A/n': f"{r['alone_A']:,}/{r['alone_n']:,}", 'Human-only runs alone: %': r['alone_percent'],
             'Human identifications (P): with / alone': f"{r['with_P']:,} / {r['alone_P']:,}"} for r in passes]
    return [('', pd.DataFrame(rows))]


# Sheets of the published workbook (D1 to D11, numbered in the order the Supplement first cites them, then D12, the
# source data of main-text Figure 2), each the tables it holds, in block order, and the number of leading label
# columns that identify a row in each block.
SHEETS = {'D1': [('s9a', 2), ('s49', 1)], 'D2': [('s61', 2)], 'D3': [('s64', 3)], 'D4': [('s62', 3), ('s63', 3)],
          'D5': [('s65', 3)], 'D6': [('s1b', 1), ('s4', 3)], 'D7': [('s51', 2)], 'D8': [('s73', 2)],
          'D9': [('s60', 3)], 'D10': [('s58', [2, 3]), ('s53', 2), ('s54', 3)], 'D11': [('s59', 6)],
          'D12': [('m3', 2), ('m4', 1)]}
# v5.1 tables cited by the short Supplement but not held in the workbook; written to a separate check workbook.
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


# ---------------------------------------------------------------------------------------------------------------
# Publication layout (carried from the build of the first published workbook; values now come from build()).

REPOSITORY = 'https://github.com/pcarvalho75/mbr-second-pass-audit'
README_TITLE = ('Supplementary Data 1. Full-detail tables behind the Supplementary Information of "Completeness at a cost: '
                'auditing the second pass of match-between-runs in single-cell DIA proteomics".')
README_NOTES = [
    'Each sheet holds one table, named by its table ID (D1 to D12); a merged table holds labelled blocks. Values are recomputed '
    'from the search outputs by the analysis code (scripts/supplementary_data.py) and written with the digits printed in the '
    'v5.1 Supplement tables and their machine-readable versions (D12: those of Tables 3 and 4 of the earlier manuscript).',
    'Plain numbers are numeric cells whose number format reproduces the printed digits; numbers with more than 15 significant '
    'digits, intervals, count pairs and labels are stored as text so that every digit is kept.',
    'Vocabulary of the older labels: second-pass-only = added; first-pass-only = lost; validation = the 250 pg dataset '
    '(PXD056832) with all 54 runs analysed together; recipient = a cell, or a run of the 250 pg dataset; on = original search '
    '(refinement on); off or no_refine = refinement off; standard or standard_library = the standard-library search; H, Y, E = '
    'human, yeast, E. coli; D, F, R, T, L = first pass, final, retained, added, lost.',
    f'Analysis code: {REPOSITORY}. Search outputs, libraries, aggregate tables and verification records: available from the corresponding author on request.',
    'The other tables of the earlier Supplement (mass accuracy, checksums, support strata, seed checks and summary tables) are '
    f'available on request, and the scripts at {REPOSITORY} regenerate them.',
]
CONFNOTE = ('Columns: original = original search (refinement on); no_refine = refinement off; standard_library = '
            'standard-library search. ')
# sheet: (title, old v5.1 numbers, Supplement section, description, note, block labels)
PUBLICATION = {
    'D1': ('Candidate ratios of the predicted and empirical libraries', 'Table S9A, Table S49', 'S2 (the candidate ratio r)',
           'Exact modified-sequence-plus-charge candidates by species class and the yeast/human ratio r, for the original '
           'predicted library and each empirical library.', '',
           ['Block A. Precursor candidate universes and charge strata of the empirical and original predicted libraries (v5.1 Table S9A)',
            'Block B. Empirical-library candidate census, refinement on and off (v5.1 Table S49)']),
    'D2': ('Accounting of the three single-cell configurations under F1, F2 and F3', 'Table S61', 'S4 (Table S4)',
           'Per-configuration counts of first-pass, final, retained, added (second_pass_only) and lost observations, both lost '
           'subclasses (library present or excluded) and the F2 reclassification counts.',
           CONFNOTE + 'Full precision, as in the machine-readable v5.1 table.', [None]),
    'D3': ('Class composition, yeast fractions, ratios to retained, conditional scales and group representation, three '
           'single-cell configurations, all filters', 'Table S64', 'S2, S4, S5 (Table S3)',
           'For every filter, class and configuration: observations, human, yeast, species-exclusive denominator, observed yeast '
           'fraction, ratio to retained, historical and empirical-library r and scales, and the share of observations whose '
           'protein group stays in the same cell (also through another peptide).',
           CONFNOTE + 'Fractions are proportions, not percentages. conditional_historical uses r = 0.276123; '
           "conditional_configuration_empirical uses the configuration's empirical-library r. Classes: first_pass, final, "
           'retained, second_pass_only (added), lost, lost_library_present, lost_library_excluded.', [None]),
    'D4': ('Class-fraction and class-ratio intervals of the three single-cell configurations (F1)', 'Table S62, Table S63',
           'S2, S4 (Table S3)', 'Point estimates, 95% limits and draw states per class or contrast, resampling scheme '
           '(recipient = cell, precursor, two_way) and configuration.',
           CONFNOTE + 'Fractions are proportions. 2,000 draws, seed 20261018.',
           ['Block A. Class-fraction intervals (v5.1 Table S62)', 'Block B. Class-ratio intervals (v5.1 Table S63)']),
    'D5': ('Represented entities and completeness of the three single-cell configurations under F1, F2 and F3', 'Table S65',
           'S2, S4 (Table S4); main-text Figure 1B',
           'Represented precursors, peptides and protein groups per pass, common rows, positive cell-entity counts, denominators, '
           'and own-row and common-row completeness. The original column under F1, F2 and F3 is the source of Figure 1B of the '
           'main text.',
           CONFNOTE + 'first_all_row_completeness and final_all_row_completeness are own-row completeness; common-row '
           'completeness uses common_row_denominator. Completeness values are proportions.', [None]),
    'D6': ('Cell-line pools of one or four lines: definitions and per-pool results', 'Table S1B, Table S4', 'S6 (Table S6)',
           'Pool sizes, arms (H homogeneous, M mixed), compositions and membership rules, with added observations and yeast per '
           'fixed pool.', 'Each row of Block B is a fixed pool, not an independent replicate.',
           ['Block A. Pool definitions (v5.1 Table S1B)', 'Block B. Per-pool results (v5.1 Table S4)']),
    'D7': ('Per-run accounting of the 18 single-organism runs, 250 pg dataset, all 54 runs (F1)', 'Table S51', 'S7 (Table S7)',
           'A, P, n, excluded assignments and observed wrong-organism fraction per run and class, including both lost subclasses.',
           'Runs 31 to 33 and 46 to 48 are human-only, 34 to 39 yeast-only and 40 to 45 E. coli-only. Fractions are proportions. '
           'Per-run intervals are unavailable.', [None]),
    'D8': ('Within-organism completeness filters, 250 pg dataset (final result, F1)', 'Table S73', 'S7 (Table S8)',
           'Rows (distinct precursors) and events (precursor-run observations) per organism with and without the 3 of 6 and 5 '
           'of 6 filters.',
           'For row units, run resampling under-represents rows seen in few runs, so run-based and two-way row intervals sit '
           'below their points; the precursor-only interval is the labelled post hoc sensitivity. Event units are unaffected.',
           [None]),
    'D9': ("Each organism's six single-organism runs analysed alone (post hoc): contrasts and class fractions", 'Table S60',
           'S8 (Table S9)', 'Block A: class-fraction and support contrasts under F1, F2 and F3; Block B: class absent fractions.',
           'Each cohort is one DIA-NN 2.5.1 two-pass search of six single-organism runs against the same three-organism library, '
           'with --no-refine-q and --report-decoys, reusing the 54-run first-pass caches. Support counts the other five runs. '
           'Intervals: two-way run and key percentile, 2,000 shared draws. D first pass, F final, R retained, T added, L lost. '
           'Source lines refer to the producer output files.',
           ['Block A. Class-fraction and support contrasts. Count pairs are A/P; fractions and their contrasts use A/n, where '
            'n=A+P. Ratios and interval limits are displayed to three decimal places; fraction differences to six decimal '
            'places. Source lines are physical TSV lines, including the header.',
            'Block B. Class absent fractions at organism scope. D is the first-pass set, F the final set, R retained, T '
            'second-pass-only and L first-pass-only. Fractions are displayed to six decimal places.']),
    'D10': ('Fold changes by context and batch, with paired eligible-key intersections (F1, Precursor.Normalised)',
            'Table S58, Table S53, Table S54', 'S10 (Table S11)',
            'Median log2 fold change, median bias, MAE, RMSE and eligible keys per context for both batches, and MAE on '
            'precursors eligible in two contexts.',
            'The batches come from separate searches and are never pooled. D_F_intersection: eligible in both first-pass and '
            'final contexts; R_T_intersection: eligible in both retained and added contexts.',
            ['Block A. Fold changes by acquisition batch (v5.1 Table S58)',
             'Block B. Paired eligible-key contexts in the 9 January replication (v5.1 Table S58, continued)',
             'Block C. 16 January contexts with universe keys, expected values and RMSE (v5.1 Table S53)',
             'Block D. Paired eligible-key intersections, 16 January (v5.1 Table S54)']),
    'D11': ('Paired, call-consequence and abundance-stratified fold-change analyses with sensitivities (post hoc)', 'Table S59',
            'S10 (Table S11)', 'P1 paired retained/final, P2 call consequences and P3 abundance terciles, primary and '
            'sensitivity rows.',
            'P1: on keys eligible in both the retained and final contexts with at least one positive added observation, the '
            'final minus retained log2 fold change, the fraction moving toward or away from zero, and the paired change in '
            "absolute error. P2: each context's own eligible keys; a missed call is an absolute log2 fold change below 1, a "
            'wrong sign the opposite sign to the design, within 0.5 inclusive. P3 (16 January only, by prior scope): abundance '
            'terciles assigned once on the union of retained and added eligible keys from the median log2 Precursor.Normalised '
            'over the twelve runs; shared keys count in both classes. D first pass, F final, R retained, T added. Intervals: '
            '2,000-draw precursor bootstraps conditional on the observed runs. Internal hash-stamped plan fixed before '
            'computation; run once; not pre-registered.',
            [None, 'Continued: P2 primary results. D denotes first pass, F final, R retained and T second-pass-only; a dash '
             'marks an inapplicable count.', 'Continued: P3 primary results.',
             'Continued: Sensitivity results, by endpoint, batch, filter (F1, F2, F3), quantity and species.']),
    'D12': ('Source data of Figure 2 (250 pg dataset, F1)', 'main-text Tables 3 and 4', 'S7, S8; main-text Figure 2',
            'Block A: identifications exclusive to an absent (A) and to the intended (P) organism, n and the observed '
            'wrong-organism fraction with its two-way 95% interval, by intended organism and class, all 54 runs analysed '
            'together (Figure 2A). Block B: the six human-only runs analysed with all 54 runs and alone, A/n and percentage '
            'for the first pass, the final result and each class, and the human identifications P of both analyses (Figure 2B).',
            'Blocks A and B were Tables 3 and 4 of the earlier manuscript and keep their printed digits. n = A + P, excluding '
            'shared and ambiguous assignments; pooled rows sum over the three organisms; intervals are two-way percentile '
            'intervals. The first pass is retained plus lost observations, and the final result retained plus added.',
            ['Block A. All 54 runs analysed together (Figure 2A; earlier Table 3)',
             'Block B. Human-only runs analysed with all 54 runs and alone, post hoc (Figure 2B; earlier Table 4)']),
}
# Display digits of the published tables, by sheet and column (every number inside the cell text is written so).
# g4: four significant digits; fN: N decimal places; thousands: integer with thousands separators; percent_or_f4:
# two decimal places for percentages, four otherwise. Columns not listed are written exactly as computed.
FORMATS = {
    'D1': {'Total keys': 'thousands', 'Human': 'thousands', 'Yeast': 'thousands', 'Shared': 'thousands', 'I/L shared': 'thousands',
           'r': 'f9', 'Yeast/human ratio': 'f12'},
    'D6': {'Second-pass-only observations': 'thousands', 'Yeast/second-pass-only observation,%': 'f3'},
    'D7': {'Observed wrong-organism fraction [two-way 95% interval]': 'g4'},
    'D8': {'Observed wrong-organism fraction [two-way 95% interval]': 'g4', 'Precursor-only 95% interval (rows; post hoc)': 'g4'},
    'D9': {'Fraction ratio': 'f3', '95% two-way interval': 'f3', 'Fraction difference': 'f6', 'absent fraction': 'f6'},
    'D10': {c: 'g4' for c in ['Jan 16 median', 'Jan 16 bias', 'Jan 16 MAE', 'Jan 9 median', 'Jan 9 bias', 'Jan 9 MAE',
                              'Median log2 fold change', 'Median bias', 'MAE', 'Expected', 'Observed median', 'RMSE']},
    'D11': {'Point': 'percent_or_f4', '95% interval': 'percent_or_f4'},
    'D12': {'Observed wrong-organism fraction [two-way 95% interval]': 'g4', 'With all 54 runs: %': 'f3', 'Human-only runs alone: %': 'f3'},
}
NUMBER_TOKEN = re.compile(r'-?\d+(?:\.\d+)?(?:e[-+]?\d+)?')


def _rule(rule, number, text):
    x = float(number)
    if rule == 'g4':
        return format(x, '.4g')
    if rule == 'thousands':
        return f'{int(x):,}'
    if rule == 'percent_or_f4':
        return format(x, '.2f' if '%' in text else '.4f')
    return format(x, '.' + rule[1:] + 'f')


def display(value, rule=None):
    """The published text of one value."""
    if isinstance(value, (bool, np.bool_)):
        return str(value)
    if isinstance(value, (int, np.integer)):
        return _rule(rule, value, '') if rule else str(int(value))
    if isinstance(value, (float, np.floating)):
        return _rule(rule, value, '') if rule else repr(float(value))
    text = str(value)
    if rule:
        text = NUMBER_TOKEN.sub(lambda m: _rule(rule, m.group(0), text), text)
    return text


INT_RE = re.compile(r'-?\d+')
COMMA_RE = re.compile(r'-?\d{1,3}(,\d{3})+')
DEC_RE = re.compile(r'-?\d+\.\d+')


def _significant(s):
    return len(s.lstrip('-').replace('.', '').replace(',', '').lstrip('0'))


def to_cell(s):
    """(value, number format) that displays the text s exactly, or (s, None) for text (as in the first build)."""
    if INT_RE.fullmatch(s) and _significant(s) <= 15 and not (s.startswith('-0') or (s.startswith('0') and len(s) > 1)):
        return int(s), '0'
    if COMMA_RE.fullmatch(s) and _significant(s) <= 15:
        return int(s.replace(',', '')), '#,##0'
    if DEC_RE.fullmatch(s) and _significant(s) <= 15:
        v = float(s)
        if v == 0 and s.startswith('-'):
            return s, None
        integer = s.lstrip('-').split('.')[0]
        if len(integer) > 1 and integer.startswith('0'):
            return s, None
        return v, '0.' + '0' * len(s.split('.')[1])
    return s, None


def shown(value, fmt):
    """The text a cell displays (read-back check)."""
    if fmt is None:
        return value
    if fmt == '0':
        return str(value)
    if fmt == '#,##0':
        return f'{value:,}'
    return f'{value:.{len(fmt) - 2}f}'


def publish(sheets, path):
    """Write the published workbook (README and D1 to D12) and check that every cell displays its text. A sheet whose
    step outputs are missing is written with the reason only (as on synthetic data)."""
    from openpyxl import Workbook, load_workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    bold, fill = Font(bold=True), PatternFill('solid', fgColor='D9D9D9')
    wb = Workbook()
    ws0 = wb.active
    ws0.title = 'README'
    ws0.append([README_TITLE])
    ws0['A1'].font = Font(bold=True, size=12)
    for note in README_NOTES:
        ws0.append([note])
    ws0.append([])
    ws0.append(['Sheet', 'Title', 'Old v5.1 number', 'Supplement section (short table)', 'Description'])
    for c in ws0[ws0.max_row]:
        c.font = bold
        c.fill = fill
    for name, (title, old, section, description, note, labels) in PUBLICATION.items():
        ws0.append([name, title, old, section, description])
    for column, width in zip('ABCDE', [8, 70, 16, 18, 90]):
        ws0.column_dimensions[column].width = width
    for row in ws0.iter_rows(min_row=1, max_row=ws0.max_row):
        for c in row:
            c.alignment = Alignment(wrap_text=True, vertical='top')
    ledger = []
    for name, (title, old, section, description, note, labels) in PUBLICATION.items():
        blocks = sheets[name]
        ws = wb.create_sheet(name)
        if isinstance(blocks, Exception):
            ws.append([f'not produced: {blocks}'])
            continue
        assert len(blocks) == len(labels), name
        ws.append([f'Table {name}. {title}'])
        ws['A1'].font = Font(bold=True, size=12)
        ws.append([f'Old v5.1 number: {old}. Supplement section: {section}. {description}'])
        if note:
            ws.append([f'Note: {note}'])
        ws.append([])
        widths, first_header = {}, None
        rules = FORMATS.get(name, {})
        for (_, frame), label in zip(blocks, labels):
            if label:
                ws.append([label])
                ws.cell(row=ws.max_row, column=1).font = Font(italic=True)
            header = list(frame.columns)
            ws.append(header)
            first_header = first_header or ws.max_row
            for c in ws[ws.max_row]:
                c.font = bold
                c.fill = fill
                c.alignment = Alignment(wrap_text=True, vertical='top')
            for record in frame.itertuples(index=False, name=None):
                texts = [display(v, rules.get(c)) for c, v in zip(header, record)]
                cells = [to_cell(s) for s in texts]
                ws.append([v for v, f in cells])
                r = ws.max_row
                for j, (s, (v, f)) in enumerate(zip(texts, cells), start=1):
                    ws.cell(row=r, column=j).number_format = f if f else '@'
                    ledger.append((name, r, j, s))
                    widths[j] = max(widths.get(j, 0), len(s))
            for j, h in enumerate(header, start=1):
                widths[j] = max(widths.get(j, 0), min(len(h), 24))
            ws.append([])
        for j, w in widths.items():
            ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 8), 60)
        if first_header:
            ws.freeze_panes = ws.cell(row=first_header + 1, column=1)
    wb.save(path)
    back = load_workbook(path)
    bad = 0
    for name, r, j, s in ledger:
        c = back[name].cell(row=r, column=j)
        fmt = None if c.number_format in ('@', 'General') and isinstance(c.value, str) else c.number_format
        if (shown(c.value, fmt) if fmt else c.value) != s:
            bad += 1
    if bad:
        raise RuntimeError(f'{bad} cells of {path} do not display their text')
    return len(ledger)


def write_extra(sheets, path):
    """The v5.1 tables behind the short Supplement tables S3, S5, S6 and S10 (check workbook, full precision)."""
    with pd.ExcelWriter(path, engine='openpyxl') as writer:
        for name in EXTRA:
            blocks = sheets[name]
            if isinstance(blocks, Exception):
                pd.DataFrame({'not produced': [str(blocks)]}).to_excel(writer, sheet_name=name[:31], index=False)
                continue
            row = 0
            for title, frame in blocks:
                frame.to_excel(writer, sheet_name=name[:31], index=False, startrow=row)
                row += len(frame) + 2
