"""Figure 2: wrong-organism identifications by observation class in the 250 pg dataset. (A) Retained, added and
lost observations of the human-only, yeast-only and E. coli-only runs, all 54 runs analysed together, with two-way
95% intervals. (B) The six human-only runs analysed with the other 48 runs and on their own: first pass, retained,
added, lost and final result.
Needs step_pg250.py --study primary and step_within_organism.py --organism H.
Outputs: <output_dir>/figures/figure2.{pdf,svg,png,tiff} and figure2-data.tsv (every plotted value with its A, P and n).
"""
from _setup import setup


def main():
    setup(__doc__)
    from mbraudit import config
    from mbraudit.figures import figure2
    root = config.output_root() / 'pg250'
    out = config.output_dir('figures')
    frame = figure2.panel_data(root / 'primary' / 'export', root / 'within-H')
    figure2.write_data(frame, out / 'figure2-data.tsv')
    figure2.draw(frame, out)
    print('Figure 2:', len(frame), 'plotted values;', out, flush=True)


if __name__ == '__main__':
    main()
