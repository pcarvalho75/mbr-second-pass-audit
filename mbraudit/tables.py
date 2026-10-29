"""Assemble the main Tables 1 and 2, and the rows behind Figure 2 (Table D12), from the pipeline outputs.

Nothing here estimates anything new: each cell is read from a pipeline output and, where the table prints a
derived figure, derived with the carried function that defined it (scale_point, final_single_cell_bound,
accepted_passes). Every function returns a list of row dictionaries with full-precision values; the entry
scripts write them as TSV and as Markdown with the manuscript's rounding. wrong_organism() and
human_only_passes() were Tables 3 and 4 of the earlier manuscript; they now feed Figure 2 (panels A and B) and
Table D12 of Supplementary Data 1.
"""
import csv
from pathlib import Path

import pandas as pd

CONFIGS = ['on', 'off', 'standard']


def _read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', keep_default_na=True)


def _rows(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


# ---------------------------------------------------------------- Table 1

def table1(single_cell):
    single_cell = Path(single_cell)
    measures = [('First-pass observations', 'acc', 'first_observations'), ('Final observations', 'acc', 'final_observations'),
                ('Retained observations', 'acc', 'retained_observations'),
                ('Added observations (final only)', 'acc', 'second_pass_only_observations'),
                ('Lost observations (first pass only)', 'acc', 'lost_observations'),
                ('Lost, in empirical library', 'acc', 'lost_library_present'),
                ('Lost, outside empirical library', 'acc', 'lost_library_excluded'),
                ('Represented precursors, first pass', 'cov', 'first_represented_entities'),
                ('Represented precursors, final', 'cov', 'final_represented_entities'),
                ('First-pass common-row completeness, %', 'cov', 'first_common_row_completeness'),
                ('Final common-row completeness, %', 'cov', 'final_common_row_completeness')]
    values = {}
    for conf in CONFIGS:
        acc = _read(single_cell / conf / 'accounting' / 'accounting.tsv').set_index('filter').loc['F1']
        cov = _read(single_cell / conf / 'accounting' / 'represented_coverage.tsv')
        cov = cov[cov['filter'].eq('F1') & cov.level.eq('precursor')].iloc[0]
        values[conf] = {'acc': acc, 'cov': cov}
    rows = []
    for label, source, column in measures:
        row = {'measure': label}
        for conf in CONFIGS:
            v = values[conf][source][column]
            row[conf] = float(v) * 100 if column.endswith('completeness') else int(v)
        rows.append(row)
    return rows


# ---------------------------------------------------------------- Table 2

def search_space_scales(candidate_counts, single_cell):
    """Restates the loop of the old analysis/search_space_scales.py main() with the carried scale_point."""
    from mbraudit.singlecell.scales import scale_point
    single_cell = Path(single_cell)
    tables = {'on': (single_cell / 'comparison-off' / 'class_species_scales_comparison.tsv', 'original_'),
              'off': (single_cell / 'comparison-off' / 'class_species_scales_comparison.tsv', 'norefine_'),
              'standard': (single_cell / 'comparison-standard' / 'class_species_scales_comparison.tsv', 'norefine_')}
    opportunities = _rows(candidate_counts)
    result = []
    for setting, (path, prefix) in tables.items():
        library = 'standard' if setting == 'standard' else 'restricted'
        selected = [r for r in opportunities if r['library'] == library and r['policy'] == 'historical_yeast_IL' and r.get('charge', 'all') == 'all']
        if len(selected) != 1:
            raise ValueError('Expected one one-way candidate ratio for ' + library)
        h, y = int(selected[0]['H']), int(selected[0]['Y'])
        for row in _rows(path):
            if row['filter'] != 'F1' or row['class'] not in ['retained', 'second_pass_only', 'lost', 'lost_library_present', 'lost_library_excluded']:
                continue
            f = row[prefix + 'observed_yeast_fraction']
            fraction, percent, shown = scale_point(f, h, y, setting, row['class'])
            result.append({'configuration': setting, 'class': row['class'], 'full_precision_f': f,
                           'human_candidates': h, 'yeast_candidates': y, 'r_exact': f'{y}/{h}',
                           'search_space_scale_fraction': fraction, 'search_space_scale_percent': percent,
                           'Search-space scale': shown,
                           'empirical_fraction_unchanged': row[prefix + 'conditional_configuration_empirical']})
    if len(result) != 15:
        raise ValueError('Expected five classes in each of three separate configurations')
    return result


def table2(single_cell, candidate_counts):
    from mbraudit.singlecell.scales import final_single_cell_bound
    single_cell = Path(single_cell)
    scales = {(r['configuration'], r['class']): r for r in search_space_scales(candidate_counts, single_cell)}
    counts = _rows(candidate_counts)
    rows, whole = [], []
    for conf in CONFIGS:
        species = _rows(single_cell / conf / 'accounting' / 'class_species.tsv')
        intervals = _read(single_cell / conf / 'intervals' / 'class_intervals.tsv')
        ratios = _read(single_cell / conf / 'intervals' / 'class_ratio_intervals.tsv')
        f1 = {r['class']: r for r in species if r['filter'] == 'F1'}
        for cls, label in [('retained', 'Retained'), ('second_pass_only', 'Added'), ('lost', 'Lost')]:
            point = intervals[intervals['class'].eq(cls) & intervals.scheme.eq('two_way')].iloc[0]
            if cls == 'retained':
                ratio = 1.0
            else:
                rr = ratios[ratios.contrast.eq(cls + '_vs_retained') & ratios.scheme.eq('two_way')].iloc[0]
                ratio = float(rr.point)
            rows.append({'class': label, 'configuration': conf, 'yeast': int(f1[cls]['yeast']),
                         'n': int(f1[cls]['species_exclusive_denominator']),
                         'observed_percent': 100 * float(f1[cls]['observed_yeast_fraction']),
                         'interval_low_percent': 100 * float(point.low95), 'interval_high_percent': 100 * float(point.high95),
                         'ratio_to_retained': ratio,
                         'upper_bound_percent': float(scales[conf, cls]['search_space_scale_percent'])})
        library = 'standard' if conf == 'standard' else 'restricted'
        item = [r for r in counts if r['library'] == library and r['policy'] == 'historical_yeast_IL' and r.get('charge', 'all') == 'all'][0]
        bound = final_single_cell_bound(species, int(item['H']), int(item['Y']))
        whole.append({'configuration': conf, **bound, 'r': int(item['Y']) / int(item['H'])})
    return rows, whole


# ---------------------------------------------------------------- Figure 2A, Table D12 Block A (earlier Table 3)

def wrong_organism(export):
    """Wrong-organism identifications of the single-organism runs, all 54 runs analysed together (F1): A, P, n,
    observed fraction A/n and two-way 95% interval per intended organism (and pooled) and class."""
    export = Path(export)
    counts = _read(export / 'endpoint_counts.tsv')
    intervals = _read(export / 'fraction_intervals.tsv')
    two_way = intervals[intervals.scheme.eq('two_way')].set_index('endpoint_id')
    rows = []
    for group, label in [('H', 'Human'), ('Y', 'Yeast'), ('E', 'E. coli'), ('all', 'Pooled')]:
        scope = 'pooled' if group == 'all' else 'organism'
        for cls, name in [('R', 'Retained'), ('T', 'Added'), ('L', 'Lost')]:
            row = counts[counts.endpoint.eq('T1') & counts['mask'].eq('F1') & counts.scope.eq(scope) & counts['group'].eq(group) & counts['class'].eq(cls)].iloc[0]
            i = two_way.loc[int(row.endpoint_id)]
            rows.append({'organism': label, 'class': name, 'A': int(row.A), 'P': int(row.P), 'n': int(row.n),
                         'fraction': float(row.absent_fraction), 'low95': float(i.low95), 'high95': float(i.high95)})
    return rows


# ---------------------------------------------------------------- Figure 2B, Table D12 Block B (earlier Table 4)

def human_only_passes(primary_export, within_human):
    """The six human-only runs analysed with the other 48 runs (all 54) and on their own (F1): A, n, percentage
    wrong-organism and P for the first pass, the final result and each class; and the added-count comparison."""
    from mbraudit.pg250.summary import accepted_passes, added_count_comparison
    mixed = _rows(Path(primary_export) / 'endpoint_counts.tsv')
    alone = _rows(Path(within_human) / 'endpoint_counts.tsv')

    def classes(rows):
        out = {}
        for cls in ['R', 'T', 'L']:
            r = [x for x in rows if x['endpoint'] == 'T1' and x['mask'] == 'F1' and x['scope'] == 'organism' and x['group'] == 'H' and x['tercile'] == 'all' and x['class'] == cls][0]
            out[cls] = (int(r['A']), int(r['P']))
        return out
    with_mix, on_own = classes(mixed), classes(alone)
    passes = {'with': {p['pass']: p for p in accepted_passes(mixed, 'H')}, 'alone': {p['pass']: p for p in accepted_passes(alone, 'H')}}
    rows = []
    for label, key in [('First pass', 'first'), ('Final result', 'final')]:
        a, b = passes['with'][key], passes['alone'][key]
        rows.append({'class': label, 'with_A': a['A'], 'with_n': a['n'], 'with_percent': 100 * a['A'] / a['n'],
                     'alone_A': b['A'], 'alone_n': b['n'], 'alone_percent': 100 * b['A'] / b['n'], 'with_P': a['P'], 'alone_P': b['P']})
    for label, cls in [('Retained', 'R'), ('Added', 'T'), ('Lost', 'L')]:
        (a1, p1), (a2, p2) = with_mix[cls], on_own[cls]
        rows.append({'class': label, 'with_A': a1, 'with_n': a1 + p1, 'with_percent': 100 * a1 / (a1 + p1),
                     'alone_A': a2, 'alone_n': a2 + p2, 'alone_percent': 100 * a2 / (a2 + p2), 'with_P': p1, 'alone_P': p2})
    comparison = added_count_comparison(mixed, alone, 'H')
    return rows, comparison
