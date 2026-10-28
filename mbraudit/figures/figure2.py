"""Figure 2: fold changes of retained and added peptide ions against their intensity (16 January mixtures).

points() is a plain function of the fold-change table and the accepted events:
- y is each eligible ion's class fold change (fc_precursors.tsv, F3, Precursor.Normalised, contexts R and T):
  median log2 quantity in 5H4Y1E minus median log2 quantity in 5H1Y4E, with at least two positive
  quantities in each mixture (computed upstream by fold_change.fc);
- x is the median log2 of the ion's positive, finite final-result Precursor.Normalised over the twelve
  5H4Y1E and 5H1Y4E runs (F3 accepted events), divided by log2(10). This is the abundance that
  posthoc_fc.analyse uses for its P3 terciles, so an ion has the same x in both panels.
Numbers are parsed with Python float() from text, as the post hoc analysis does, so values agree to the last bit.
draw() renders the two panels with 0.3-log10 binned medians (bins with at least 10 ions).
"""
import math
from pathlib import Path

import numpy as np
import pandas as pd

QUANTITY = 'Precursor.Normalised'
COLORS = {'Y': '#E69F00', 'E': '#CC79A7'}
TITLES = {'R': 'retained', 'T': 'added'}
PANELS = {'R': 'A', 'T': 'B'}
BIN_WIDTH = .3
MIN_BIN_N = 10


def _text(path, usecols=None):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, usecols=usecols)


def points(fc_path, keys_path, events_path, samples_path):
    fc = _text(fc_path)
    fc = fc[(fc['mask'] == 'F3') & (fc.quantity == QUANTITY) & fc.context.isin(['R', 'T']) & (fc.eligible == 'True')].copy()
    fc['key'] = fc.key.map(float).astype(int)
    fc['observed_log2fc'] = fc.observed_log2fc.map(float)
    keys = _text(keys_path)
    keys['key'] = keys.key.map(float).astype(int)
    samples = _text(samples_path)
    axis = samples[(samples.kind == 'mixture') & samples['group'].isin(['5H4Y1E', '5H1Y4E'])]
    events = _text(events_path, usecols=['sample', 'Modified.Sequence', 'Precursor.Charge', 'species', 'mask', 'stage', QUANTITY])
    events = events[(events.stage == 'final') & (events['mask'] == 'F3') & events['sample'].isin(axis['sample'])].copy()
    events[QUANTITY] = events[QUANTITY].map(lambda x: np.nan if x == 'NA' else float(x))
    events = events[events[QUANTITY].gt(0) & np.isfinite(events[QUANTITY])]
    events = events.merge(keys[['key', 'Modified.Sequence', 'Precursor.Charge']], on=['Modified.Sequence', 'Precursor.Charge'], how='left', validate='many_to_one')
    events['log_abundance'] = np.log2(events[QUANTITY])
    abundance = events.groupby(['species', 'key']).log_abundance.median().rename('abundance_log2').reset_index()
    frame = fc.merge(abundance, on=['species', 'key'], how='left', validate='many_to_one')
    if frame.abundance_log2.isna().any():
        raise ValueError('An eligible ion has no positive final quantity in the twelve runs')
    frame['log10_intensity'] = frame.abundance_log2 / np.log2(10.)
    frame = frame.merge(keys[['key', 'Modified.Sequence', 'Precursor.Charge']], on='key', how='left', validate='many_to_one')
    columns = ['species', 'context', 'key', 'Modified.Sequence', 'Precursor.Charge', 'n_high', 'n_low', 'observed_log2fc', 'abundance_log2', 'log10_intensity']
    return frame[columns].sort_values(['context', 'species', 'key']).reset_index(drop=True)


def binned_medians(x, y, edges):
    bins = np.digitize(x, edges) - 1
    rows = []
    for i in range(len(edges) - 1):
        take = bins == i
        n = int(take.sum())
        valid = n >= MIN_BIN_N
        rows.append({'left': float(edges[i]), 'right': float(edges[i + 1]), 'n': n, 'shown': valid,
                     'x_median': float(np.median(x[take])) if valid else float((edges[i] + edges[i + 1]) / 2),
                     'y_median': float(np.median(y[take])) if valid else np.nan})
    return rows


def draw(frame, output, name='figure2'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    from PIL import Image
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.labelsize': 8, 'axes.titlesize': 9,
                         'xtick.labelsize': 7.2, 'ytick.labelsize': 7.2, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .6, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'svg.hashsalt': 'fold-change-dotplot'})
    fig, axes = plt.subplots(1, 2, figsize=(174 / 25.4, 80 / 25.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.095, right=.99, bottom=.19, top=.86, wspace=.12)
    xlo, xhi = frame.log10_intensity.min() - .08, frame.log10_intensity.max() + .08
    ylo = math.floor((min(frame.observed_log2fc.min(), -2) - .2) * 2) / 2
    yhi = math.ceil((max(frame.observed_log2fc.max(), 2) + .2) * 2) / 2
    edges = np.arange(math.floor(xlo / BIN_WIDTH) * BIN_WIDTH, math.ceil(xhi / BIN_WIDTH) * BIN_WIDTH + BIN_WIDTH / 2, BIN_WIDTH)
    trends = []
    for ax, ctx in zip(axes, ['R', 'T']):
        for yy, color in [(2, COLORS['Y']), (0, '#A8ADB3'), (-2, COLORS['E'])]:
            ax.axhline(yy, color=color, lw=.7, ls=(0, (4, 3)), alpha=.85, zorder=0)
        for sp in ['Y', 'E']:
            p = frame[(frame.context == ctx) & (frame.species == sp)]
            x = p.log10_intensity.to_numpy(float)
            y = p.observed_log2fc.to_numpy(float)
            ax.scatter(x, y, s=5, alpha=.36, color=COLORS[sp], edgecolors='none', rasterized=False, zorder=2)
            rows = binned_medians(x, y, edges)
            trends.extend({'context': ctx, 'species': sp, **r} for r in rows)
            ax.plot([r['x_median'] for r in rows], [r['y_median'] for r in rows], color=COLORS[sp], lw=1.05, zorder=3)
        ax.set_title(TITLES[ctx], pad=8)
        ax.text(-.02, 1.035, PANELS[ctx], transform=ax.transAxes, ha='right', va='bottom', fontsize=10, fontweight='bold')
        if ctx == 'R':
            ax.text(xlo + .06, yhi - .75, 'yeast', color=COLORS['Y'], fontsize=8, fontweight='bold', ha='left', va='center')
            ax.text(xlo + .06, ylo + .65, 'E. coli', color=COLORS['E'], fontsize=8, fontweight='bold', fontstyle='italic', ha='left', va='center')
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.tick_params(width=.6, length=3, pad=3)
    fig.supxlabel('log10 intensity', fontsize=8, y=.035)
    fig.supylabel('log2 fold change', fontsize=8, x=.017)
    fig.savefig(output / (name + '.pdf'), metadata={'Creator': None, 'Producer': None, 'CreationDate': None, 'ModDate': None})
    fig.savefig(output / (name + '.svg'), metadata={'Date': None, 'Creator': None})
    fig.savefig(output / (name + '.png'), dpi=600)
    plt.close(fig)
    with Image.open(output / (name + '.png')) as image:
        image.convert('RGB').save(output / (name + '.tiff'), compression='tiff_lzw', dpi=(600, 600))
    pd.DataFrame(trends).to_csv(output / (name + '-binned-medians.tsv'), sep='\t', index=False, lineterminator='\n', na_rep='NA')
    return trends
