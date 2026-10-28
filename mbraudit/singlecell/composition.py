"""Pools of one cell line against pools of four lines: yeast per added observation.

Carried from analysis/revision_identification.py (composition and helpers), analysis/composition.py
(rate_ratio_ci). OUT is set by the entry script (plumbing).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from mbraudit.species import clopper_pearson
from mbraudit.singlecell.reports import Reports, KEY, H6_RUNS, M6, H16, M16, natural, transfers, count_record, SEED, DRAWS
from mbraudit.singlecell.bootstrap import bootstrap, interval
OUT = None


SCHEMES = ['recipient', 'precursor', 'two_way']


def checkpoint(done, total):
    if done % 500 == 0:
        print(f'bootstrap {done}/{total}', flush=True)


def save(name, frame):
    pd.DataFrame(frame).to_csv(OUT / name, sep='\t', index=False)


def uncertainty(name, frames, cells, kind, point):
    rows = []
    for scheme in SCHEMES:
        print(name, scheme, flush=True)
        values = bootstrap(frames, cells, kind, scheme, checkpoint=checkpoint)
        save(f'draws_{name}_{scheme}.tsv', {'draw': np.arange(1, DRAWS + 1), 'value': values})
        row = dict(cohort=name, scheme=scheme, point=point, **interval(values))
        row.update(recipients=len(cells), precursor_blocks=len(set().union(*(set(f.key) for f in frames))))
        rows.append(row)
    return rows


def per_cell(frame, cells):
    return pd.DataFrame([dict(cell=cell, **count_record(frame[frame.Run.eq(cell)])) for cell in cells])


def composition(reports):
    intervals = []
    cells_out = []
    pools = []
    paired = []
    arms = {}
    for n, h_runs, m_runs in [(6, list(H6_RUNS), list(M6)), (16, list(H16), list(M16))]:
        cells = sorted(set((c for members in (M6 if n == 6 else M16).values() for c in members)), key=natural)
        h = transfers(reports.arm(h_runs))
        m = transfers(reports.arm(m_runs))
        arms[n, 'H'] = h
        arms[n, 'M'] = m
        hr = count_record(h)
        mr = count_record(m)
        point, lo, hi = rate_ratio_ci(mr['yeast'], mr['events'], hr['yeast'], hr['events'])
        intervals.append({'cohort': f'composition_{n}', 'scheme': 'conditional_count_model', 'point': point, 'low95': lo, 'high95': hi, 'draws': 0, 'undefined_draws': 0, 'infinite_draws': 0, 'interval_status': 'conditional count model', 'recipients': len(cells), 'precursor_blocks': len(set(h.key) | set(m.key))})
        intervals.extend(uncertainty(f'composition_{n}', [h, m], cells, 'composition', point))
        per = {}
        for label, frame in [('H', h), ('M', m)]:
            pc = per_cell(frame, cells)
            pc['pool_size'] = n
            pc['composition'] = label
            pc['observed_yeast_per_transfer'] = np.divide(pc.yeast, pc.events, out=np.full(len(pc), np.nan), where=pc.events.gt(0))
            cells_out.extend(pc.to_dict('records'))
            per[label] = pc.set_index('cell')
            for pool, g in frame.groupby('pool'):
                record = count_record(g)
                pools.append(dict(pool_size=n, composition=label, pool=pool, recipients=g.Run.nunique(), observed_yeast_per_transfer=record['yeast'] / record['events'], **record))
        for cell in cells:
            hh, mm = (per['H'].loc[cell], per['M'].loc[cell])
            paired.append({'pool_size': n, 'cell': cell, 'homogeneous_events': int(hh.events), 'mixed_events': int(mm.events), 'homogeneous_yeast': int(hh.yeast), 'mixed_yeast': int(mm.yeast), 'homogeneous_yeast_per_transfer': hh.observed_yeast_per_transfer, 'mixed_yeast_per_transfer': mm.observed_yeast_per_transfer, 'paired_yeast_rate_difference': mm.observed_yeast_per_transfer - hh.observed_yeast_per_transfer, 'homogeneous_estimated_combined_fdp': hh.estimated_combined_fdp, 'mixed_estimated_combined_fdp': mm.estimated_combined_fdp, 'paired_combined_fdp_difference': mm.estimated_combined_fdp - hh.estimated_combined_fdp})
    save('A_composition_uncertainty.tsv', intervals)
    save('A_per_cell.tsv', cells_out)
    save('A_per_pool.tsv', pools)
    save('A_paired_differences.tsv', paired)
    return (intervals, arms)


def _theta_from_pi(pi: float, n1: float, n2: float) -> float:
    if pi <= 0:
        return 0.0
    if pi >= 1:
        return float('inf')
    return pi / (1.0 - pi) * (n2 / n1)


def rate_ratio_ci(e1: int, n1: int, e2: int, n2: int, alpha: float=0.05) -> tuple[float, float, float]:
    e1, n1, e2, n2 = (int(e1), int(n1), int(e2), int(n2))
    if n1 <= 0 or n2 <= 0 or (e1 == 0 and e2 == 0):
        return (float('nan'), float('nan'), float('nan'))
    if e2 == 0:
        point = float('inf')
    else:
        point = e1 / n1 / (e2 / n2)
    plo, phi = clopper_pearson(e1, e1 + e2)
    return (point, _theta_from_pi(plo, n1, n2), _theta_from_pi(phi, n1, n2))
