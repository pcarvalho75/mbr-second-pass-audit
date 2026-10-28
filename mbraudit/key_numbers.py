"""Every number quoted in the manuscript text, read from the pipeline outputs.

compute(root) returns rows of (id, value, source). Ids are listed with their manuscript location and printed
value in references/published-values.tsv; scripts/verify.py compares the two.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from mbraudit import tables


def _read(path, **kw):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', **kw)


def _one(frame, **criteria):
    sel = frame
    for k, v in criteria.items():
        sel = sel[sel[k].astype(str) == str(v)]
    if len(sel) != 1:
        raise ValueError(f'Expected one row for {criteria}, found {len(sel)}')
    return sel.iloc[0]


def compute(root, candidate_counts):
    root = Path(root)
    sc, pg = root / 'single-cell', root / 'pg250'
    out = []

    def add(key, value, source):
        out.append({'id': key, 'value': value, 'source': source})

    # Table 1 and the single-cell accounting numbers.
    for row in tables.table1(sc):
        for conf in tables.CONFIGS:
            add(f"T1|{row['measure']}|{conf}", row[conf], f'single-cell/{conf}/accounting')
    on = sc / 'on' / 'accounting'
    cls = _read(on / 'class_species.tsv')
    lost = _one(cls, filter='F1', **{'class': 'lost'})
    add('P44|lost kept group via other peptide', int(lost.same_cell_group_other_peptide_in_final), 'class_species.tsv F1 lost')
    add('P44|lost kept group via other peptide, %', 100 * float(lost.same_cell_group_other_peptide_in_final_fraction), 'class_species.tsv F1 lost')
    ent = _read(on / 'entity_species.tsv')
    add('P44|precursors wholly absent', int(_one(ent, filter='F1', level='precursor', state='wholly_absent').observations), 'entity_species.tsv')
    hum = _read(on / 'human_group_losses.tsv')
    add('P44|human-only protein groups wholly absent', int(_one(hum, filter='F1').human_groups_wholly_absent), 'human_group_losses.tsv')
    cov = _read(on / 'represented_coverage.tsv')
    for filt in ['F1', 'F2', 'F3']:
        row = _one(cov, filter=filt, level='precursor')
        for phase in ['first', 'final']:
            for universe in ['all', 'common']:
                add(f'Fig1B|{filt}|{phase}|{universe}', 100 * float(row[f'{phase}_{universe}_row_completeness']), 'represented_coverage.tsv')
    # Table 2 and the trap numbers.
    rows, whole = tables.table2(sc, candidate_counts)
    for r in rows:
        for field in ['yeast', 'n', 'observed_percent', 'interval_low_percent', 'interval_high_percent', 'ratio_to_retained', 'upper_bound_percent']:
            add(f"T2|{r['class']}|{r['configuration']}|{field}", r[field], f"single-cell/{r['configuration']}")
    for w in whole:
        add(f"T2|final|{w['configuration']}|observed_percent", 100 * float(w['lower_fraction']), 'final_single_cell_bound')
        add(f"T2|final|{w['configuration']}|upper_bound_percent", 100 * float(w['estimated_upper_fraction']), 'final_single_cell_bound')
        add(f"T2|final|{w['configuration']}|yeast/n", f"{w['yeast']}/{w['denominator']}", 'final_single_cell_bound')
        add(f"T2|r|{w['configuration']}", w['r'], 'candidate counts, historical_yeast_IL')
    ratios = _read(sc / 'on' / 'intervals' / 'class_ratio_intervals.tsv')
    for name in ['second_pass_only_vs_retained', 'lost_vs_retained']:
        row = _one(ratios, contrast=name, scheme='two_way')
        add(f'P49|{name}|low95', float(row.low95), 'class_ratio_intervals.tsv')
        add(f'P49|{name}|high95', float(row.high95), 'class_ratio_intervals.tsv')
    pools = _read(sc / 'pools' / 'A_composition_uncertainty.tsv')
    for n in [16, 6]:
        add(f'P51|mixed/single-line ratio, pools of {n}', float(_one(pools, cohort=f'composition_{n}', scheme='two_way').point), 'A_composition_uncertainty.tsv')
        add(f'S6|two-way interval, pools of {n}', f"{float(_one(pools, cohort=f'composition_{n}', scheme='two_way').low95):.2f} to {float(_one(pools, cohort=f'composition_{n}', scheme='two_way').high95):.2f}", 'A_composition_uncertainty.tsv')
    sup = _read(sc / 'on' / 'intervals' / 'support_contrasts.tsv')
    row = _one(sup, cohort='all_cells_64', filter='F1', bin_scheme='current_three_positive', quantity_tercile='high', support_bin='low_positive')
    add('P70|single-cell support ratio', float(row.fraction_ratio), 'support_contrasts.tsv')
    add('P70|single-cell support counts', f'{int(row.numerator_yeast)}/{int(row.numerator_n)} vs {int(row.reference_yeast)}/{int(row.reference_n)}', 'support_contrasts.tsv')

    # The 250 pg dataset.
    export = pg / 'primary' / 'export'
    for r in tables.table3(export):
        for field in ['A', 'P', 'n', 'fraction', 'low95', 'high95']:
            add(f"T3|{r['organism']}|{r['class']}|{field}", r[field], 'pg250/primary/export')
    counts = _read(export / 'endpoint_counts.tsv', keep_default_na=False)
    design = pd.read_csv(Path(__file__).resolve().parents[1] / 'design' / 'samples54.tsv', sep='\t', keep_default_na=False)
    human_runs = design[(design.kind == 'single') & (design.present_species == 'H')]['sample']
    t1 = counts[(counts.endpoint == 'T1') & (counts['mask'] == 'F1') & (counts.scope == 'run')]
    added = t1[(t1['class'] == 'T') & t1['group'].isin(human_runs)].set_index('group').absent_fraction.astype(float)
    kept = t1[(t1['class'] == 'R') & t1['group'].isin(human_runs)].set_index('group').absent_fraction.astype(float)
    add('P53|per-run added fraction min, %', 100 * added.min(), 'endpoint_counts.tsv run scope')
    add('P53|per-run added fraction max, %', 100 * added.max(), 'endpoint_counts.tsv run scope')
    add('P53|added above retained in every human-only run', bool((added > kept.reindex(added.index)).all()), 'endpoint_counts.tsv run scope')
    t4, comparison = tables.table4(export, pg / 'within-H')
    for r in t4:
        for field in ['with_A', 'with_n', 'with_percent', 'alone_A', 'alone_n', 'alone_percent', 'with_P', 'alone_P']:
            add(f"T4|{r['class']}|{field}", r[field], 'accepted_passes / endpoint_counts')
    add('P56|fewer wrong-organism additions alone, %', float(comparison['wrong_added_count_reduction_percent']), 'added_count_comparison')
    from mbraudit.pg250.summary import accepted_passes
    rows = [dict(r) for r in counts.astype(str).to_dict('records')]
    for organism in ['Y', 'E']:
        first = [p for p in accepted_passes(rows, organism) if p['pass'] == 'first'][0]
        add(f'P54|first-pass wrong-organism fraction {organism}, %', 100 * first['A'] / first['n'], 'accepted_passes')
    points = _read(export / 'contrasts.tsv')
    limits = _read(export / 'contrast_intervals.tsv')
    for numerator in ['T', 'L']:
        sel = dict(endpoint='T1', mask='F1', scope='pooled', group='all', numerator=numerator, reference='R')
        p = points.loc[(points[list(sel)].astype(str) == pd.Series(sel)).all(axis=1)].iloc[0]
        lim = limits.loc[(limits[list(sel)].astype(str) == pd.Series(sel)).all(axis=1) & limits.scheme.eq('two_way') & limits.statistic.eq('ratio')].iloc[0]
        add(f'P54|pooled {numerator}/R ratio', float(p.fraction_ratio), 'contrasts.tsv')
        add(f'P54|pooled {numerator}/R low95', float(lim.low95), 'contrast_intervals.tsv')
        add(f'P54|pooled {numerator}/R high95', float(lim.high95), 'contrast_intervals.tsv')
    sel = dict(endpoint='T2', mask='F1', scope='pooled', group='all', tercile='high', numerator='positive_60_100', reference='positive_0_20')
    p = points.loc[(points[list(sel)].astype(str) == pd.Series(sel)).all(axis=1)].iloc[0]
    add('P70|250 pg support ratio high/low', float(p.fraction_ratio), 'contrasts.tsv T2')
    add('P70|250 pg support counts', f'{int(p.numerator_A)}/{int(p.numerator_n)} vs {int(p.reference_A)}/{int(p.reference_n)}', 'contrasts.tsv T2')
    # Completeness filter (Table S8, old S73).
    fold = pg / 'primary' / 'fold-change'
    t3 = _read(fold / 't3_summary.tsv', keep_default_na=False)
    cov3 = _read(fold / 't3_coverage.tsv', keep_default_na=False)
    row = _one(t3, mask='F1', organism='H', scheme='stage_specific', stage='final', threshold=5, transition='F', scope='group', unit='events')
    add('P56|human final wrong-organism fraction, at least 5 of 6, %', 100 * float(row.false_fraction), 't3_summary.tsv')
    row = _one(cov3, mask='F1', organism='H', scheme='stage_specific', stage='final', threshold=5, transition='F')
    add('P56|human identifications kept, at least 5 of 6, %', 100 * float(row.present_event_coverage), 't3_coverage.tsv')
    # Fold changes (Results 3.5, Figure 2).
    post = _read(pg / 'posthoc-fc' / 'summary.tsv', keep_default_na=False)
    base = post[(post.batch == 'jan16') & (post['mask'] == 'F1') & (post.quantity == 'Precursor.Normalised')]
    for sp in ['Y', 'E']:
        for ctx in ['R', 'T']:
            row = _one(base, species=sp, endpoint='P2', context=ctx, metric='missed_twofold')
            add(f'P61|{sp}|{ctx}|missed two-fold', f'{int(float(row.numerator))}/{int(float(row.eligible_keys))}', 'posthoc summary P2')
            add(f'P61|{sp}|{ctx}|missed two-fold, %', 100 * float(row.point), 'posthoc summary P2')
        for ctx in ['R', 'T']:
            row = _one(base, species=sp, endpoint='P2', context=ctx, metric='wrong_sign')
            add(f'P61|{sp}|{ctx}|wrong sign, %', 100 * float(row.point), 'posthoc summary P2')
        low = base[(base.species == sp) & (base.endpoint == 'P3') & (base.tercile == 'low')]
        add(f'P62|{sp}|low-tercile T minus R bias', float(_one(low, metric='median_bias_T_minus_R').point), 'posthoc summary P3')
        add(f'P62|{sp}|added ions in low tercile', int(float(_one(low, metric='median_bias_T_minus_R').T_keys)), 'posthoc summary P3')
        add(f'P62|{sp}|P1 MAE change', float(_one(base, species=sp, endpoint='P1', metric='MAE_change').point), 'posthoc summary P1')
    for study, tag in [('primary', 'jan16'), ('replication', 'jan09')]:
        fs = _read(pg / study / 'fold-change' / 'fc_summary.tsv', keep_default_na=False)
        fs = fs[(fs['mask'] == 'F1') & (fs.quantity == 'Precursor.Normalised')]
        for ctx in ['D', 'R', 'F', 'T']:
            add(f'P61|{tag}|yeast median log2 FC|{ctx}', float(_one(fs, species='Y', context=ctx, population='own_context').median_observed_log2fc), 'fc_summary.tsv')
        better = 0
        for sp in ['Y', 'E']:
            r = _one(fs, species=sp, context='R', population='R_T_intersection')
            t = _one(fs, species=sp, context='T', population='R_T_intersection')
            better += int(float(t.MAE) < float(r.MAE))
            add(f'P62|{tag}|{sp}|R_T intersection keys', int(r.eligible_keys), 'fc_summary.tsv')
        add(f'P62|{tag}|added MAE smaller on R/T intersection (of 2)', better, 'fc_summary.tsv')
    fig2 = root / 'figures' / 'figure2-points.tsv'
    if fig2.exists():
        f2 = _read(fig2)
        for sp in ['Y', 'E']:
            for ctx in ['R', 'T']:
                add(f'Fig2|{sp}|{ctx}|ions', int(((f2.species == sp) & (f2.context == ctx)).sum()), 'figure2-points.tsv')
            shared = set(f2[(f2.species == sp) & (f2.context == 'R')].key) & set(f2[(f2.species == sp) & (f2.context == 'T')].key)
            add(f'Fig2|{sp}|ions in both panels', len(shared), 'figure2-points.tsv')
    qa = pg / 'native-qa' / 'native-run-qa.tsv'
    if qa.exists():
        q = _read(qa, keep_default_na=False)
        for _, r in q[q.kind == 'blank'].iterrows():
            add(f"S10|{r['pass']}|{r.manifest_ordinal}|identified precursors", int(r['Precursors.Identified']), 'native-run-qa.tsv')
    carry = pg / 'carryover' / 'carryover.tsv'
    if carry.exists():
        c = _read(carry, keep_default_na=False, dtype={'acquisition_id': str, 'organism': str})
        for r in c.itertuples():
            add('|'.join(['S9', r.section] + [x for x in [r.acquisition_id, r.organism] if x] + [r.measure]), r.value, 'carryover.tsv')
    add('P40|designed runs', len(design), 'design/samples54.tsv')
    for kind, n in design.kind.value_counts().items():
        add(f'P40|runs of kind {kind}', int(n), 'design/samples54.tsv')
    return out
