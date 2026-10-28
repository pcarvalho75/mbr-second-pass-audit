"""Single-cell dataset: side-by-side tables of the original search (refinement on) and each repeat search
(refinement off; standard library), and each configuration's empirical-library candidate ratio (Tables 1, 2, S72).

Needs step_single_cell.py for on, off and standard. Outputs: <output_dir>/single-cell/comparison-off and
comparison-standard (prefix original_ = refinement on, norefine_ = the repeat search).
"""
import os

from _setup import setup, note


def main():
    setup(__doc__)
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    from mbraudit import config
    from mbraudit.singlecell import configurations as cfg
    from mbraudit.singlecell.reports import Reports

    base = config.output_root() / 'single-cell'
    census = base / 'empirical-census-on.tsv'
    reports = Reports()
    for repeat in ['off', 'standard']:
        os.environ['MBRAUDIT_SINGLE_CELL'] = repeat
        comparison = config.output_dir('single-cell', 'comparison-' + repeat)
        oldr, newr = cfg.empirical_opportunities(reports, census, comparison)
        cfg.comparison_tables(base / 'on' / 'accounting', base / 'on' / 'intervals',
                              base / repeat / 'accounting', base / repeat / 'intervals', comparison, oldr, newr)
        note(comparison / 'step.json', {'step': 'compare-configurations', 'original': 'on', 'repeat': repeat,
                                        'empirical_r_original': oldr, 'empirical_r_repeat': newr,
                                        'repeat_library': str(config.single_cell_report_dir('all_cells_64') / 'report-lib.parquet')})
        print('Comparison on versus', repeat, 'complete:', comparison, flush=True)


if __name__ == '__main__':
    main()
