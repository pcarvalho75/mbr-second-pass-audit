"""Figure 1 panel B source and shared drawing helpers (carried from figures/make_figures.py).

TABLES and OUTPUT are set by the entry script (plumbing).
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
TABLES = None
OUTPUT = None


BLUE, ORANGE, GREY, BLACK = ('#0072B2', '#D55E00', '#777777', '#222222')


plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 13, 'axes.labelsize': 13, 'axes.titlesize': 14, 'xtick.labelsize': 13, 'ytick.labelsize': 13, 'legend.fontsize': 12.5, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none', 'svg.hashsalt': 'single-cell-dia-revision'})


INPUTS = {}


VALUES = []


FIGURE_ORDER = {1: 'S5'}


FIGURE_NAMES = {'figure01_composition': 'figureS5_composition'}


def read(name, json_data=False):
    path = TABLES / name
    raw = path.read_bytes()
    INPUTS[name] = hashlib.sha256(raw).hexdigest()
    if json_data:
        return json.loads(raw)
    return list(csv.DictReader(raw.decode('utf-8-sig').splitlines(), delimiter='\t'))


def one(rows, **criteria):
    selected = [r for r in rows if all((r[k] == str(v) for k, v in criteria.items()))]
    if len(selected) != 1:
        raise ValueError(f'Expected one row for {criteria},found {len(selected)}')
    return selected[0]


def record(figure, label, row, fields):
    VALUES.append({'figure': FIGURE_ORDER.get(figure, figure), 'label': label, 'values': {k: row[k] for k in fields}})


def value(row, key):
    return float(row[key]) if row[key] not in ['', None] else math.nan


def interval(ax, y, row, point, low, high, color=BLUE, scale=1, marker='o'):
    x, lo, hi = [value(row, k) * scale for k in [point, low, high]]
    ax.plot([lo, hi], [y, y], color=color, linewidth=1.8)
    ax.plot([lo, lo], [y - 0.07, y + 0.07], color=color)
    ax.plot([hi, hi], [y - 0.07, y + 0.07], color=color)
    ax.plot(x, y, marker=marker, color=color, markersize=6, linestyle='none')


def save(fig, name):
    name = FIGURE_NAMES.get(name, name)
    fig.tight_layout(pad=1.2, h_pad=2, w_pad=2)
    png = OUTPUT / f'{name}.png'
    fig.savefig(png, dpi=300, bbox_inches='tight', metadata={'Author': '[AUTHOR]', 'Title': name})
    fig.savefig(OUTPUT / f'{name}.pdf', bbox_inches='tight', metadata={'Author': '[AUTHOR]', 'Creator': '[AUTHOR]', 'CreationDate': None, 'ModDate': None})
    fig.savefig(OUTPUT / f'{name}.svg', bbox_inches='tight', metadata={'Date': None, 'Creator': '[AUTHOR]'})
    svg = OUTPUT / f'{name}.svg'
    svg.write_bytes(svg.read_bytes().replace(b'\r\n', b'\n'))
    plt.close(fig)
    with Image.open(png) as image:
        image.convert('RGB').save(OUTPUT / f'{name}.tiff', compression='tiff_lzw', dpi=(300, 300))


def primary_flow():
    rows = read('r4a_figure_flow.tsv')
    accounting = read('r4a_central_accounting.tsv')
    coverage = read('r4a_figure_completeness.tsv')
    base = one(accounting, filter='F1', level='precursor')
    selected = {name: one(rows, filter='F1', **{'class': name}) for name in ['first_pass', 'retained', 'second_pass_only', 'lost_library_present', 'lost_library_excluded']}
    fig = plt.figure(figsize=(11.5, 10.2))
    grid = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.3, 1.9], width_ratios=[2.2, 1.0], hspace=0.4, wspace=0.55)
    a = fig.add_subplot(grid[0, :])
    details = fig.add_subplot(grid[1, 0])
    loss = fig.add_subplot(grid[1, 1])
    b = fig.add_subplot(grid[2, :])
    colors = {'retained': BLUE, 'second_pass_only': '#009E73', 'lost_library_present': ORANGE, 'lost_library_excluded': GREY}
    labels = {'retained': 'Retained', 'second_pass_only': 'Second-pass-only', 'lost_library_present': 'First-pass-only:\nin library', 'lost_library_excluded': 'First-pass-only:\noutside library'}
    for y, names in [(1, ['retained', 'lost_library_present', 'lost_library_excluded']), (0.53, ['retained', 'second_pass_only'])]:
        start = 0
        for name in names:
            count = int(selected[name]['observations'])
            a.barh(y, count, left=start, height=0.17, color=colors[name])
            start += count
    maximum = int(base['final_observations']) * 1.06
    class_rows = []
    for name, row in selected.items():
        if name == 'first_pass':
            continue
        record('1', name, row, ['observations', 'yeast', 'species_exclusive_denominator', 'observed_yeast_pct', 'same_cell_group_other_peptide_pct', 'group_annotation_evaluable_observations'])
        class_rows.append([labels[name].replace('First-pass-only:\n', 'Loss: ').replace('\n', ' '), f"{int(row['observations']):,}", f"{float(row['observed_yeast_pct']):.3f}", f"{float(row['same_cell_group_other_peptide_pct']):.2f}"])
    details.axis('off')
    table = details.table(cellText=class_rows, colLabels=['Observation class', 'Count', 'Yeast\n(%)', 'Group*\n(%)'], cellLoc='center', colWidths=[0.43, 0.2, 0.16, 0.21], loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(14.5)
    table.scale(1, 2.3)
    for (i, j), cell in table.get_celld().items():
        cell.set_edgecolor('#D9D9D9')
        cell.set_linewidth(0.7)
        if i == 0:
            cell.set_facecolor('#ECECEC')
            cell.set_text_props(weight='bold')
    for i, name in enumerate(['retained', 'second_pass_only', 'lost_library_present', 'lost_library_excluded'], 1):
        table[i, 0].get_text().set_color(colors[name])
    loss_labels = ['In library', 'Outside library']
    for y, name in enumerate(['lost_library_present', 'lost_library_excluded']):
        count = int(selected[name]['observations'])
        loss.barh(y, count, color=colors[name], height=0.43)
        loss.text(count + 600, y, f'{count:,}', va='center', fontsize=13)
        loss.text(0, y - 0.3, loss_labels[y], va='bottom', fontsize=13)
    loss.set(yticks=[], ylim=(1.65, -0.65), xlim=(0, 41000), xlabel='First-pass-only observations', title='Expanded loss classes')
    loss.tick_params(labelsize=12)
    for key in ['first_observations', 'final_observations']:
        record('1', key, base, [key])
    a.set(yticks=[1, 0.53], yticklabels=[f"First pass\n{int(base['first_observations']):,}", f"Final\n{int(base['final_observations']):,}"], ylim=(0.2, 1.35), xlim=(0, maximum), title='A  Primary F1 accepted-set accounting', xlabel='Precursor-cell observations')
    from matplotlib.ticker import FuncFormatter
    loss.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{v / 1000:g}k' if v else '0'))
    a.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{v / 1000:g}k' if v else '0'))
    for j, filt in enumerate(['F1', 'F2', 'F3']):
        row = one(coverage, filter=filt, level='precursor')
        for delta, phase, universe, color, hatch in [(-0.36, 'first', 'all', BLUE, '//'), (-0.12, 'final', 'all', ORANGE, '//'), (0.12, 'first', 'common', BLUE, ''), (0.36, 'final', 'common', ORANGE, '')]:
            field = phase + '_' + universe + '_row_completeness_pct'
            v = value(row, field)
            label = ('First pass' if phase == 'first' else 'Final') + (' own rows' if universe == 'all' else ' common rows')
            b.barh(j + delta, v, height=0.15, color=color, hatch=hatch, edgecolor='white', linewidth=0.7, label=label if j == 0 else None)
            b.text(v + 0.3, j + delta, f'{v:.2f}', va='center', fontsize=14)
            record('1', filt + ' ' + field, row, [field])
    b.set(yticks=range(3), yticklabels=['F1 (primary)', 'F2', 'F3 (naive)'], ylim=(2.6, -0.6), xlim=(0, 56), xlabel='Finite positive precursor-cell entries (%)', title='B  Own-row and common-row completeness')
    b.legend(frameon=False, ncol=2, loc='upper left', bbox_to_anchor=(0, -0.22), fontsize=13.5)
    save(fig, 'figure01_accounting')
