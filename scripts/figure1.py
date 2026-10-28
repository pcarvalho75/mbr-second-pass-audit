"""Figure 1: (A) accepted observations split into retained, added and lost, with represented precursors, for the
three single-cell configurations; (B) own-row and common-row completeness under F1, F2 and F3; (C) class / retained
ratios of the wrong-organism fraction with two-way 95% intervals, single cells and the 250 pg dataset.
Needs step_single_cell.py (all three), step_pg250.py --study primary.
Outputs: <output_dir>/figures/figure1.{svg,pdf,png,tiff} and figure1-data.json (every plotted value).
"""
from _setup import setup


def main():
    setup(__doc__)
    from mbraudit import config
    from mbraudit.figures.figure1 import make
    root = config.output_root()
    outputs = make(root / 'single-cell', root / 'pg250' / 'primary' / 'export', config.output_dir('figures'))
    print('Figure 1:', ', '.join(sorted(outputs)), flush=True)


if __name__ == '__main__':
    main()
