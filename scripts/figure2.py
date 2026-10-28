"""Figure 2: fold changes of retained (A) and added (B) yeast and E. coli peptide ions between the 5H4Y1E and
5H1Y4E mixtures (16 January, six runs each) against their intensity, with 0.3-log10 binned medians.
Needs step_pg250.py --study primary. Outputs: <output_dir>/figures/figure2.* and figure2-points.tsv (every dot).
"""
from _setup import setup


def main():
    setup(__doc__)
    from mbraudit import config
    from mbraudit.figures import figure2
    root = config.output_root() / 'pg250' / 'primary'
    out = config.output_dir('figures')
    frame = figure2.points(root / 'fold-change' / 'fc_precursors.tsv', root / 'fold-change' / 'keys.tsv',
                           root / 'export' / 'accepted_events.tsv', config.resolve(config.value('pg250.primary.design')))
    frame.to_csv(out / 'figure2-points.tsv', sep='\t', index=False, float_format='%.17g', lineterminator='\n')
    figure2.draw(frame, out)
    print('Figure 2:', len(frame), 'ions;', out, flush=True)


if __name__ == '__main__':
    main()
