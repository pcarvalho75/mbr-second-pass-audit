"""Carryover or false matches (post hoc; Supplement S9): the wrong-organism identifications added to the six
human-only runs of the 54-run 250 pg search, against the yeast-only and E. coli-only runs acquired earlier the
same day.

Uses the package's own definitions: report rows are targets (Decoy = 0) named by endpoints.canonical_sample,
accepted under filters.filter_mask with the F1 fields of endpoints.MASKS (missing or non-finite values fail), keyed
by endpoints.KEY, called by organism from the fixed species map of the three-organism predicted library (as
endpoints.adapt_reports does), and classed by endpoints.memberships (R retained, T added). Acquisition order comes
from the acquisition_id column of the design. Four checks, all descriptive:
- run order: added yeast and E. coli identifications per human-only run, in acquisition order;
- organism: E. coli against yeast additions;
- intensity: each added ion's Precursor.Quantity over its median Precursor.Quantity in the same-day source runs of
  its organism (final result, F1), the median and interquartile range of the log10 ratios, and the Pearson
  correlation of the log10 intensities;
- recurrence: in how many of the six runs each added wrong-organism precursor is added.
Also the accepted target precursors of the blank runs, by organism, in each pass.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from mbraudit.filters import filter_mask
from mbraudit.pg250.endpoints import EVENT, KEY, MASKS, canonical_sample, memberships

MASK = 'F1'
SPECIES_EXCLUSIVE = ['H', 'Y', 'E']
COLUMNS = ['Run', 'Modified.Sequence', 'Precursor.Charge', 'Decoy', 'Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'Precursor.Quantity']


def accepted(report, design, species_map, stage, mask=MASK):
    """Accepted target observations of one report under one filter, with their organism call."""
    raw = pd.read_parquet(report, columns=COLUMNS)
    raw['sample'] = raw.Run.map(canonical_sample)
    if not set(raw['sample']) <= set(design['sample']):
        raise ValueError('Unexpected report runs')
    raw['Precursor.Charge'] = pd.to_numeric(raw['Precursor.Charge'], errors='raise').astype(int)
    raw = raw[raw.Decoy.eq(0)].drop(columns=['Run']).drop_duplicates()
    if raw.duplicated(EVENT).any():
        raise ValueError('Conflicting duplicate run/key rows')
    q = raw[MASKS[mask]].apply(pd.to_numeric, errors='coerce')
    keep = filter_mask(q.fillna(2.), mask) & np.isfinite(q).all(axis=1)
    out = raw.loc[keep, ['sample'] + KEY + ['Precursor.Quantity']].copy()
    out['Precursor.Quantity'] = pd.to_numeric(out['Precursor.Quantity'], errors='coerce')
    out = out.merge(species_map, on=KEY, how='left', validate='many_to_one')
    out['species'] = out.species.fillna('unknown')
    out['mask'], out['stage'] = mask, stage
    return out


def analyse(first, final, design, native_qa=None, mask=MASK):
    """Per-run and summary values (one long table) from the accepted first-pass and final observations."""
    events = pd.concat([first, final], ignore_index=True)
    status = pd.DataFrame([{'sample': s, 'stage': st, 'mask': mask, 'status': 'complete'} for s in design['sample'] for st in ['first', 'final']])
    members, _ = memberships(events, status, design, None)
    members = members[members['mask'].eq(mask)]
    acq = dict(zip(design['sample'], design.acquisition_id.astype(int)))
    singles = design[design.kind.eq('single')]
    human = singles[singles.present_species.eq('H')].copy()
    human['acq'] = human['sample'].map(acq)
    human = human.sort_values('acq')
    day = set(human.filename_date.astype(str))
    rows = []

    def add(section, measure, value, sample='', organism=''):
        rows.append({'section': section, 'sample': sample, 'acquisition_id': acq.get(sample, ''), 'organism': organism, 'measure': measure, 'value': value})

    exclusive = members[members['sample'].isin(human['sample']) & members.species.isin(SPECIES_EXCLUSIVE)]
    for cls, label in [('R', 'retained'), ('T', 'added')]:
        part = exclusive[exclusive['class'].eq(cls)]
        add('check', f'{label} wrong-organism', int(part.species.ne('H').sum()))
        add('check', f'{label} species-exclusive', len(part))
    added = exclusive[exclusive['class'].eq('T')]
    for sample in human['sample']:
        g = added[added['sample'].eq(sample)]
        d = exclusive[exclusive['sample'].eq(sample) & exclusive['class'].eq('D')]
        add('run order', 'added species-exclusive', len(g), sample)
        add('run order', 'added yeast', int(g.species.eq('Y').sum()), sample, 'Y')
        add('run order', 'added E. coli', int(g.species.eq('E').sum()), sample, 'E')
        add('run order', 'added wrong-organism fraction, %', 100 * float(g.species.ne('H').mean()), sample)
        add('run order', 'first-pass wrong-organism', int(d.species.ne('H').sum()), sample)
        add('run order', 'first-pass species-exclusive', len(d), sample)
    add('organism', 'added yeast', int(added.species.eq('Y').sum()), organism='Y')
    add('organism', 'added E. coli', int(added.species.eq('E').sum()), organism='E')
    source_final = final[final['sample'].isin(singles['sample'])]
    for sp, name in [('Y', 'yeast'), ('E', 'E. coli')]:
        sources = singles[singles.present_species.eq(sp) & singles.filename_date.astype(str).isin(day)]['sample']
        med = source_final[source_final['sample'].isin(sources)].groupby(KEY)['Precursor.Quantity'].median()
        pct = med.rank(pct=True)
        a = added[added.species.eq(sp)].copy()
        k = list(zip(a['Modified.Sequence'], a['Precursor.Charge']))
        a['source_median'] = [med.get(x, np.nan) for x in k]
        a['source_percentile'] = [pct.get(x, np.nan) for x in k]
        a['ratio'] = a['Precursor.Quantity'] / a['source_median']
        seen = a.source_median.notna()
        lr = np.log10(a.ratio.dropna())
        add('intensity', f'added {name} ions', len(a), organism=sp)
        add('intensity', f'added {name} ions accepted in the source runs', int(seen.sum()), organism=sp)
        add('intensity', 'median ratio to the source median', float(10 ** lr.median()), organism=sp)
        add('intensity', 'interquartile range of log10 ratios', float(lr.quantile(.75) - lr.quantile(.25)), organism=sp)
        add('intensity', 'Pearson r of log10 intensities', float(np.corrcoef(np.log10(a.loc[seen, 'Precursor.Quantity']), np.log10(a.loc[seen, 'source_median']))[0, 1]), organism=sp)
        add('intensity', 'median source-intensity percentile', float(a.source_percentile.median()), organism=sp)
        for sample, g in a.groupby('sample', sort=False):
            r = g.ratio.dropna()
            if len(r):
                add('intensity', 'median ratio to the source median, this run', float(10 ** np.log10(r).median()), sample, sp)
    wrong = added[added.species.ne('H')].groupby(KEY).size()
    add('recurrence', 'distinct wrong-organism precursors', len(wrong))
    for n, count in wrong.value_counts().sort_index().items():
        add('recurrence', f'added in {int(n)} of the six runs', int(count))
    add('recurrence', 'added in three or more of the six runs', int(wrong.ge(3).sum()))
    blanks = design[design.kind.eq('blank')].copy()
    blanks['acq'] = blanks['sample'].map(acq)
    qa = None
    if native_qa is not None and Path(native_qa).exists():
        qa = pd.read_csv(native_qa, sep='\t', dtype=str, keep_default_na=False)
    for sample in blanks.sort_values('acq')['sample']:
        for stage, frame in [('first', first), ('final', final)]:
            g = frame[frame['sample'].eq(sample)]
            add('blanks', f'{stage} accepted target precursors', len(g), sample)
            for sp in ['H', 'Y', 'E', 'shared', 'shared_IL', 'unknown']:
                add('blanks', f'{stage} accepted {sp}', int(g.species.eq(sp).sum()), sample, sp)
            if qa is not None:
                row = qa[qa['sample'].eq(sample) & qa['pass'].eq(stage)]
                if len(row) == 1:
                    add('blanks', f'{stage} native identified precursors', int(row['Precursors.Identified'].iloc[0]), sample)
                    add('blanks', f'{stage} calibration failures', int(row.calibration_failures.iloc[0]), sample)
    return pd.DataFrame(rows, columns=['section', 'sample', 'acquisition_id', 'organism', 'measure', 'value'])


def run(first_report, final_report, predicted_species, design_path, out, native_qa=None):
    design = pd.read_csv(design_path, sep='\t', keep_default_na=False)
    species_map = pd.read_parquet(predicted_species, columns=KEY + ['species'])
    if species_map.duplicated(KEY).any():
        raise ValueError('Predicted map must have unique keys')
    first = accepted(first_report, design, species_map, 'first')
    final = accepted(final_report, design, species_map, 'final')
    table = analyse(first, final, design, native_qa)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / 'carryover.tsv', sep='\t', index=False, lineterminator='\n', float_format='%.17g')
    return table
