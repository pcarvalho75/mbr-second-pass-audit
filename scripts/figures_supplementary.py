"""Supplementary Figures S1 to S3: (S1) wrong-organism fraction of the 250 pg added identifications by first-pass
support and quantity tercile; (S2) entities wholly absent from the final result and their yeast fractions (single
cells, F1); (S3) wrong-organism fraction and intended-species coverage after the completeness filters (250 pg).
Needs step_single_cell.py --configuration on, step_pg250.py --study primary.
Outputs: <output_dir>/figures/supplementary/FigureS1.{pdf,png,tiff} to FigureS3, and figure-s-data.json (every
plotted value).
"""
from _setup import setup


def main():
    setup(__doc__)
    from mbraudit import config
    from mbraudit.figures.supplementary import make
    root = config.output_root()
    written = make(root / 'single-cell', root / 'pg250' / 'primary', config.output_dir('figures', 'supplementary'))
    print('Supplementary figures:', ', '.join(written), flush=True)


if __name__ == '__main__':
    main()
