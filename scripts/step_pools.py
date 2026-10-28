"""Single-cell dataset: pools of one cell line against pools of all four lines (Results 3.2; Supplement S6,
Table S6; Supplementary Data D6). The ratio is yeast observations per added observation, mixed over homogeneous pools, at pool
sizes 6 and 16, with a count-model interval and recipient, precursor-block and two-way bootstrap intervals.
Outputs: <output_dir>/single-cell/pools.
"""
from _setup import setup, note


LINE = {'Hela': 'HeLa', 'Hep': 'HepG2', 'A549': 'A549', 'H1975': 'H1975'}


def pool_membership(reports, comp):
    """Cells of each pool search, as named in its reports (Supplementary Data D6, Block A)."""
    import pandas as pd
    rows = []
    for alias in list(comp.H6_RUNS) + list(comp.M6) + list(comp.H16) + list(comp.M16):
        cells = sorted(set(reports.load(alias)['first'].Run) | set(reports.load(alias)['final'].Run), key=comp.natural)
        rows.extend({'pool': alias, 'cell': cell, 'line': LINE[cell.split('_')[0]]} for cell in cells)
    return pd.DataFrame(rows)


def main():
    setup(__doc__)
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    from mbraudit import config
    from mbraudit.singlecell import composition as comp
    from mbraudit.singlecell.reports import Reports

    comp.OUT = config.output_dir('single-cell', 'pools')
    reports = Reports()
    comp.composition(reports)
    pool_membership(reports, comp).to_csv(comp.OUT / 'pool_membership.tsv', sep='\t', index=False, lineterminator='\n')
    note(comp.OUT / 'step.json', {'step': 'pools', 'seed': comp.SEED, 'draws': comp.DRAWS, 'input_sha256': reports.pins})
    print('Pool composition complete:', comp.OUT, flush=True)


if __name__ == '__main__':
    main()
