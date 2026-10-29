"""Figure 2: wrong-organism identifications by observation class in the 250 pg dataset.

panel_data() collects every plotted value from the pipeline outputs through the functions that define them
(tables.wrong_organism and tables.human_only_passes); nothing is estimated here.
- Panel A, all 54 runs analysed together (F1): the observed wrong-organism fraction A/n of retained, added and lost
  observations in the human-only, yeast-only and E. coli-only runs, with the two-way 95% percentile interval
  (endpoint_counts.tsv and fraction_intervals.tsv of the 54-run analysis, T1, organism scope).
- Panel B, the six human-only runs: the same fraction for the first pass (retained plus lost), the retained, added
  and lost classes and the final result (retained plus added), analysed with the other 48 runs and on their own
  (endpoint_counts.tsv of the 54-run and of the within-human analyses).
Values are plotted as percentages (100 times the fraction). The full rows, with the pooled rows of panel A, are
Table D12 of Supplementary Data 1 (Tables 3 and 4 of the earlier manuscript).
draw() renders the figure at two-column width (7.2 in). It uses the font of Figure 1 and the type sizes Figure 1
has when printed at the same width (Figure 1 is drawn 11.945 in wide, so each of its sizes is scaled by
7.2 / 11.945), and the class colours of Figure 1.
"""
from pathlib import Path

import pandas as pd

from mbraudit import tables

COLORS = {'Retained': '#0072B2', 'Added': '#009E73', 'Lost': '#D55E00'}  # the class colours of Figure 1
INK, LINK, REFERENCE = '#333333', '#BBBBBB', '#777777'
ORGANISMS = {'Human': 'Human-only runs', 'Yeast': 'Yeast-only runs', 'E. coli': r'$\mathit{E.\:coli}$-only runs'}
CLASSES_A = ['Retained', 'Added', 'Lost']
CLASSES_B = ['First pass', 'Retained', 'Added', 'Lost', 'Final result']
LABELS_B = {'First pass': 'First\npass', 'Final result': 'Final\nresult'}  # two lines, so that labels do not touch
ANALYSES_B = {'with': 'With the other 48 runs', 'alone': 'On their own'}
TICKS_A = [0.2, 0.5, 1, 2, 5, 10, 20, 50]
TICKS_B = [0.3, 0.5, 1, 2, 5]
WIDTH, HEIGHT = 7.2, 3.2  # inches
SCALE = WIDTH / 11.945  # Figure 1 page width, inches
# Figure 1 type sizes in points: panel titles 15, axis labels and value ticks 13, category labels 12, legends 11.5;
# its interval lines are 1.8 points wide. Markers are 4.5 points across.
TITLE, LABEL, TICK, CATEGORY, LEGEND = (15 * SCALE, 13 * SCALE, 13 * SCALE, 12 * SCALE, 11.5 * SCALE)
MARKER, LINE = 4.5, 1.8 * SCALE


def panel_data(primary_export, within_human):
    """Every plotted value of Figure 2 as rows (panel, x label, class, analysis, A, P, n, percent and interval)."""
    rows = []
    for r in tables.wrong_organism(primary_export):
        if r['organism'] not in ORGANISMS:
            continue  # pooled rows are in Table D12 only
        rows.append({'panel': 'A', 'x': r['organism'] + '-only runs', 'class': r['class'], 'analysis': 'all 54 runs',
                     'A': r['A'], 'P': r['P'], 'n': r['n'], 'percent': 100 * r['fraction'],
                     'low95_percent': 100 * r['low95'], 'high95_percent': 100 * r['high95']})
    passes, _ = tables.human_only_passes(primary_export, within_human)
    passes = {r['class']: r for r in passes}
    for cls in CLASSES_B:
        r = passes[cls]
        for key, analysis in ANALYSES_B.items():
            rows.append({'panel': 'B', 'x': cls, 'class': cls, 'analysis': analysis.lower(), 'A': r[key + '_A'],
                         'P': r[key + '_P'], 'n': r[key + '_n'], 'percent': r[key + '_percent'],
                         'low95_percent': float('nan'), 'high95_percent': float('nan')})
    frame = pd.DataFrame(rows)
    if frame.duplicated(['panel', 'x', 'class', 'analysis']).any() or len(frame) != 19:
        raise ValueError('Figure 2 needs 9 panel A points and 10 panel B points')
    return frame


def write_data(frame, path):
    """The plotted values as text with every digit (repr of the float)."""
    out = frame.copy()
    for column in ['percent', 'low95_percent', 'high95_percent']:
        out[column] = [repr(float(v)) if v == v else 'NA' for v in out[column]]
    out.to_csv(path, sep='\t', index=False, lineterminator='\n')


def _percent(value, _position=None):
    return f'{value:g}%'


def _log_axis(ax, ticks):
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    ax.set_yscale('log')
    ax.set_ylim(ticks[0], ticks[-1])
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FuncFormatter(_percent))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.axhline(1, color=REFERENCE, linestyle='--', linewidth=.8, zorder=0, gid='one-percent-line')
    ax.set_ylabel('Wrong-organism identifications')


def draw(frame, output, name='figure2'):
    """Render Figure 2 (PDF, SVG, 300 dpi PNG and LZW TIFF) from panel_data()."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.transforms import blended_transform_factory
    from PIL import Image
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    style = {'font.family': 'DejaVu Sans', 'font.size': TICK, 'axes.titlesize': TITLE, 'axes.labelsize': LABEL,
             'xtick.labelsize': CATEGORY, 'ytick.labelsize': TICK, 'legend.fontsize': LEGEND, 'axes.linewidth': .6,
             'xtick.major.width': .6, 'ytick.major.width': .6, 'xtick.major.size': 3, 'ytick.major.size': 3,
             'axes.spines.top': False, 'axes.spines.right': False, 'mathtext.fontset': 'dejavusans',
             'svg.fonttype': 'none', 'pdf.fonttype': 42, 'svg.hashsalt': 'figure2-250pg'}
    with plt.rc_context(style):
        fig, (a, b) = plt.subplots(1, 2, figsize=(WIDTH, HEIGHT), gridspec_kw={'width_ratios': [1.25, 1]})
        fig.subplots_adjust(left=.085, right=.99, bottom=.135, top=.91, wspace=.3)
        panel = frame[frame.panel == 'A']
        for i, organism in enumerate(ORGANISMS):
            for j, cls in enumerate(CLASSES_A):
                r = panel[(panel.x == organism + '-only runs') & (panel['class'] == cls)].iloc[0]
                y, lo, hi = float(r.percent), float(r.low95_percent), float(r.high95_percent)
                bars = a.errorbar(i + (j - 1) * .22, y, yerr=[[y - lo], [hi - y]], fmt='o', color=COLORS[cls], markersize=MARKER,
                                  elinewidth=LINE, capsize=2.2, capthick=LINE, label=cls if i == 0 else None)
                tag = f"panelA-{organism.replace('. ', '')}-{cls}"  # SVG ids, to read the drawn values back
                bars.lines[0].set_gid(tag + '-point')
                for cap, end in zip(bars.lines[1], ['low', 'high']):
                    cap.set_gid(f'{tag}-{end}')
                bars.lines[2][0].set_gid(tag + '-interval')
        _log_axis(a, TICKS_A)
        a.text(1, 1, '1%', transform=blended_transform_factory(a.transAxes, a.transData), ha='right', va='bottom',
               color=REFERENCE, fontsize=CATEGORY)
        a.set_xlim(-.5, len(ORGANISMS) - .5)
        a.set_xticks(range(len(ORGANISMS)), list(ORGANISMS.values()))
        a.set_title('A  All 54 runs analysed together', loc='left')
        a.legend(frameon=False, loc='upper left', handletextpad=.3, borderaxespad=.2)
        panel = frame[frame.panel == 'B']
        xs = range(len(CLASSES_B))
        values = {key: [float(panel[(panel.x == cls) & (panel.analysis == analysis.lower())].iloc[0].percent) for cls in CLASSES_B]
                  for key, analysis in ANALYSES_B.items()}
        for x, w, o in zip(xs, values['with'], values['alone']):
            b.plot([x, x], [w, o], color=LINK, linewidth=2 * LINE * .65, solid_capstyle='butt', zorder=1)
        size = (MARKER * 1.1) ** 2
        b.scatter(xs, values['with'], s=size, color=INK, linewidths=.8, zorder=2, label=ANALYSES_B['with'], gid='panelB-with')
        b.scatter(xs, values['alone'], s=size, facecolor='white', edgecolor=INK, linewidths=.8, zorder=2,
                  label=ANALYSES_B['alone'], gid='panelB-alone')
        _log_axis(b, TICKS_B)
        b.set_xlim(-.5, len(CLASSES_B) - .5)
        b.set_xticks(list(xs), [LABELS_B.get(cls, cls) for cls in CLASSES_B])
        b.set_title('B  Human-only runs', loc='left')
        b.legend(frameon=False, loc='upper left', handletextpad=.3, borderaxespad=.2)
        metadata = {'Creator': None, 'Producer': None, 'CreationDate': None, 'ModDate': None}
        fig.savefig(output / (name + '.pdf'), metadata=metadata)
        fig.savefig(output / (name + '.svg'), metadata={'Date': None, 'Creator': None})
        fig.savefig(output / (name + '.png'), dpi=300)
        plt.close(fig)
    svg = output / (name + '.svg')
    svg.write_bytes(svg.read_bytes().replace(b'\r\n', b'\n'))
    with Image.open(output / (name + '.png')) as image:
        image.convert('RGB').save(output / (name + '.tiff'), compression='tiff_lzw', dpi=(300, 300))
    return [output / (name + suffix) for suffix in ['.pdf', '.svg', '.png', '.tiff']]
