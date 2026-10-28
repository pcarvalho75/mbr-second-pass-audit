"""Single-cell dataset, one search configuration: retained, added and lost observations, represented entities,
completeness, yeast trap fractions with two-way intervals, and first-pass support (Tables 1, 2, S71, S72; Figure 1).

Order of calls is that of the original analysis: accounting, support strata, workflows, candidate scales, figure
tables, then joint class intervals and support intervals.
Outputs: <output_dir>/single-cell/<configuration>/accounting and .../intervals.
"""
import os

from _setup import setup, note


def main():
    args = setup(__doc__, lambda p: p.add_argument('--configuration', choices=['on', 'off', 'standard'], required=True))
    os.environ['MBRAUDIT_SINGLE_CELL'] = args.configuration
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    from mbraudit import config
    from mbraudit.singlecell import accounting as acc, intervals, support
    from mbraudit.singlecell.census import empirical_census
    from mbraudit.singlecell.reports import Reports

    root = config.output_dir('single-cell', args.configuration)
    out4 = config.output_dir('single-cell', args.configuration, 'accounting')
    out5 = config.output_dir('single-cell', args.configuration, 'intervals')
    reports = Reports()
    # This configuration's own empirical-library census (Supplementary Data D5; the empirical scale of Table S3).
    # The original run passed the refinement-on census to every configuration; each now uses its own.
    census = config.output_root() / 'single-cell' / ('empirical-census-' + args.configuration + '.tsv')
    library = config.single_cell_report_dir('all_cells_64') / 'report-lib.parquet'
    empirical_census(reports.registry, library).to_csv(census, sep='\t', index=False, lineterminator='\n')
    reports.pin(census, 'candidate_ratio_source')
    cache = {'all_cells_64': acc.load_quantity(reports, 'all_cells_64')}
    acc.all_cell_accounting(reports, cache['all_cells_64'], out4)
    acc.support_strata(reports, cache, out4)
    acc.workflows(reports, cache, out4)
    acc.class_ratios_and_scales(out4, census)
    acc.figures(out4)
    intervals.class_intervals(reports, cache['all_cells_64'], out5)
    support.support_intervals(reports, cache['all_cells_64'], out5)
    note(root / 'step.json', {'step': 'single-cell', 'configuration': args.configuration,
                              'report_folder': str(config.single_cell_report_dir('all_cells_64')),
                              'seed': acc.SEED, 'draws': acc.DRAWS, 'input_sha256': reports.pins})
    print('Single-cell configuration', args.configuration, 'complete:', root, flush=True)


if __name__ == '__main__':
    main()
