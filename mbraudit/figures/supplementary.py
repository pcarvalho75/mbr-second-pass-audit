"""Supplementary Figures S1 to S3, from the pipeline outputs (renderers carried from the earlier figure scripts).

- Figure S1 (quantity terciles of the 250 pg added identifications) and Figure S3 (completeness selection in the
  250 pg dataset) are single panels of the earlier four-panel validation figure: render_validation is carried from
  figures/make_validation.py:render_figure, compose_validation from figures/compose_main.py (panel placement and
  letters), and split_panel from the step that cut each panel out without its letter for the short Supplement.
- Figure S2 (entities wholly absent from the final result, single cells) is carried from
  figures/make_figures.py:primary_entities, with the primary filter labelled F1.

validation_data() collects every plotted value as the text of the output tables (no statistic is computed here).
TABLES-style plumbing (panels.TABLES, panels.OUTPUT) is set by make().
"""
import json
import shutil
from pathlib import Path

import pandas as pd

# Single-panel figures of the short Supplement (name, old letter on the composed page, page clip, titles).
SPLITS = [
    {'name': 'FigureS1', 'old_panel': 'A', 'pdf_clip': [180, 0, 540, 330], 'old_title': 'A  Quantity terciles', 'title': 'Quantity terciles'},
    {'name': 'FigureS3', 'old_panel': 'B', 'pdf_clip': [0, 330, 390, 648], 'old_title': 'B  Completeness selection', 'title': 'Completeness selection'},
]
WIDTH = 720
DPI = 300


def _text(path):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)


def _one(frame, **where):
    selected = frame
    for k, v in where.items():
        selected = selected[selected[k] == str(v)]
    if len(selected) != 1:
        raise ValueError(f'Expected one row for {where}, found {len(selected)}')
    return selected.iloc[0]


def validation_data(export, fold_change, mask='F1'):
    """Plotted cells of the validation panels, copied as text from the 54-run outputs (primary filter)."""
    export, fold_change = Path(export), Path(fold_change)
    counts = _text(export / 'endpoint_counts.tsv')
    limits = _text(export / 'fraction_intervals.tsv')
    t3 = _text(fold_change / 't3_summary.tsv')
    coverage = _text(fold_change / 't3_coverage.tsv')
    fc = _text(fold_change / 'fc_summary.tsv')
    data = {'A': [], 'B': [], 'C': [], 'D': []}
    for org in ['H', 'Y', 'E', 'all']:
        for cls in ['R', 'T', 'L']:
            row = _one(counts, endpoint='T1', mask=mask, scope='pooled' if org == 'all' else 'organism', group=org, **{'class': cls})
            ci = _one(limits, endpoint_id=row['endpoint_id'], scheme='two_way')
            data['A'].append({**{k: row[k] for k in ['A', 'P', 'n', 'absent_fraction']},
                              **{k: ci[k] for k in ['low95', 'high95', 'interval_status']}, 'organism': org, 'context': cls})
    for org in ['H', 'Y', 'E']:
        for threshold in ['0', '3', '5']:
            where = dict(mask=mask, organism=org, scheme='stage_specific', stage='final', threshold=threshold, transition='F')
            row = _one(t3, **where, scope='group', unit='events')
            cov = _one(coverage, **where)
            data['C'].append({**{k: row[k] for k in ['false_fraction', 'two_way_ci_low', 'two_way_ci_high', 'A', 'P', 'n']},
                              **{k: cov[k] for k in ['present_event_coverage', 'present_events_retained', 'present_events_baseline']},
                              'organism': org, 'threshold': threshold})
    for q in ['low', 'middle', 'high', 'missing_quantity']:
        for bin_name in ['zero', 'positive_0_20', 'positive_20_40', 'positive_40_60', 'positive_60_100']:
            row = _one(counts, endpoint='T2', mask=mask, scope='pooled', group='all', tercile=q, support_bin=bin_name)
            ci = _one(limits, endpoint_id=row['endpoint_id'], scheme='two_way')
            data['B'].append({**{k: row[k] for k in ['A', 'P', 'n', 'absent_fraction']},
                              **{k: ci[k] for k in ['low95', 'high95', 'interval_status']}, 'tercile': q, 'support_bin': bin_name})
    for species in ['Y', 'E']:
        for context in ['D', 'F', 'R', 'T']:
            row = _one(fc, mask=mask, species=species, quantity='Precursor.Normalised', context=context, population='own_context')
            data['D'].append({**{k: row[k] for k in ['eligible_keys', 'universe_keys', 'expected_log2fc', 'median_observed_log2fc', 'median_bias', 'MAE', 'RMSE']},
                              'species': species, 'context': context})
    return data


CLASSES = {'R': 'Retained', 'T': 'Second-pass-only', 'L': 'First-pass-only'}
ORG = {'H': 'Human', 'Y': 'Yeast', 'E': 'E. coli'}


def num(value):
    import math
    return float(value) if value not in ['', None, 'NA', 'nan'] else math.nan


def render_validation(data, stem):
    """The four-panel validation figure (carried from figures/make_validation.py:render_figure, measured data)."""
    import math
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), layout='constrained')
    colors = {'R': '#0072B2', 'T': '#D55E00', 'L': '#009E73', 'H': '#0072B2', 'Y': '#D55E00', 'E': '#009E73'}

    def point(ax, x, row, value, lo, hi, color, missing_level=0):
        v = num(row[value])
        low = num(row[lo])
        high = num(row[hi])
        if not math.isfinite(v):
            ax.text(x, 0.025, 'x', transform=ax.get_xaxis_transform(), ha='center', fontsize=7, color=color)
            return
        has_ci = math.isfinite(low) and math.isfinite(high)
        ax.plot(x, v, 'o', color=color, ms=4, markerfacecolor=color if has_ci else 'white')
        if math.isfinite(low) and math.isfinite(high):
            ax.vlines(x, low, high, color=color, lw=1.2)
            ax.plot([x - 0.045, x + 0.045], [low, low], color=color)
            ax.plot([x - 0.045, x + 0.045], [high, high], color=color)
    ax = axes[0, 0]
    for i, org in enumerate(['H', 'Y', 'E', 'all']):
        for j, cls in enumerate(['R', 'T', 'L']):
            r = next((r for r in data['A'] if r['organism'] == org and r['context'] == cls))
            point(ax, i + (j - 1) * 0.2, r, 'absent_fraction', 'low95', 'high95', colors[cls])
    ax.set(xticks=range(4), xticklabels=['Human', 'Yeast', 'E. coli', 'Pooled'], ylabel='Observed wrong-organism fraction', title='A  Observation classes')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylim(bottom=0)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=colors[c], label=CLASSES[c]) for c in ['R', 'T', 'L']], frameon=False, fontsize=7, loc='upper right')
    ax = axes[0, 1]
    bins = ['zero', 'positive_0_20', 'positive_20_40', 'positive_40_60', 'positive_60_100']
    for j, q in enumerate(['low', 'middle', 'high']):
        for i, bin_name in enumerate(bins):
            r = next((r for r in data['B'] if r['tercile'] == q and r['support_bin'] == bin_name))
            point(ax, i + (j - 1) * 0.16, r, 'absent_fraction', 'low95', 'high95', ['#0072B2', '#009E73', '#D55E00'][j], j)
    ax.set(xticks=range(5), xticklabels=['0', '(0,.2]', '(.2,.4]', '(.4,.6]', '(.6,1]'], xlabel='First-pass support share', ylabel='Observed wrong-organism fraction', title='B  Quantity terciles', xlim=(-0.4, 4.4))
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylim(bottom=0)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=c, label=q) for q, c in zip(['Low', 'Middle', 'High'], ['#0072B2', '#009E73', '#D55E00'])], frameon=False, fontsize=7)
    ax.text(0.02, 0.98, 'x = not estimable (no plotted value)', transform=ax.transAxes, fontsize=7, va='top')
    ax = axes[1, 0]
    right = ax.twinx()
    right.spines['top'].set_visible(False)
    for j, org in enumerate(['H', 'Y', 'E']):
        coverage = []
        for i, t in enumerate(['0', '3', '5']):
            r = next((r for r in data['C'] if r['organism'] == org and r['threshold'] == t))
            point(ax, i + (j - 1) * 0.12, r, 'false_fraction', 'two_way_ci_low', 'two_way_ci_high', colors[org])
            coverage.append(num(r['present_event_coverage']))
        right.plot(range(3), coverage, '--', color=colors[org], lw=1, label=ORG[org])
    ax.set(xticks=range(3), xticklabels=['None', '3/6', '5/6'], ylabel='Observed wrong-organism fraction (events)', title='C  Completeness selection')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylim(bottom=0)
    right.set(ylabel='Intended-species event coverage', ylim=(0, 1.05))
    right.yaxis.set_major_formatter(PercentFormatter(1))
    right.legend(frameon=False, fontsize=7, loc='lower left')
    ax.text(0.03, 0.98, 'Open points: interval unavailable', transform=ax.transAxes, fontsize=7, va='top')
    ax = axes[1, 1]
    for j, sp in enumerate(['Y', 'E']):
        expected = num(next((r for r in data['D'] if r['species'] == sp))['expected_log2fc'])
        ax.axhline(expected, color=colors[sp], ls=':', lw=1)
        for i, c in enumerate(['D', 'F', 'R', 'T']):
            r = next((r for r in data['D'] if r['species'] == sp and r['context'] == c))
            x = i + (j - 0.5) * 0.15
            v = num(r['median_observed_log2fc'])
            ax.plot(x, v, 'o', color=colors[sp], ms=4)
            ax.annotate('n=' + r['eligible_keys'], (x, v), xytext=(0, 7 if j == 0 else -13), textcoords='offset points', ha='center', fontsize=6)
    ax.set(xticks=range(4), xticklabels=['First', 'Final', 'Retained', 'New'], ylabel='Median normalized log2 fold change', title='D  Mixture fold changes')
    ax.margins(y=0.14)
    ax.legend(handles=[Line2D([], [], marker='o', color=colors[s], ls=':', label=ORG[s] + '; dotted = expected') for s in ['Y', 'E']], frameon=False, fontsize=7)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix('.png'), dpi=300, metadata={'Author': '[AUTHOR]'})
    fig.savefig(stem.with_suffix('.svg'), metadata={'Creator': '[AUTHOR]', 'Date': None})
    svg = stem.with_suffix('.svg')
    svg.write_bytes(svg.read_text(encoding='utf-8').encode('utf-8'))
    fig.savefig(stem.with_suffix('.pdf'), metadata={'Author': '[AUTHOR]', 'Creator': '[AUTHOR]', 'CreationDate': None, 'ModDate': None})
    plt.close(fig)


def compose_validation(output):
    """Panel placement and letters of the composed validation page (carried from figures/compose_main.py)."""
    import fitz
    from mbraudit.figures import compose as c
    old_v = fitz.open(output / 'figure04-validation.pdf')
    doc = fitz.open()
    page = doc.new_page(width=720, height=648)
    clip = fitz.Rect(360, 0, 720, 330)
    target = fitz.Rect(180, 0, 540, 330)
    page.show_pdf_page(target, old_v, 0, clip=clip)
    c.label_pdf(page, old_v[0], clip, target, 'B  Quantity', 'A')
    clip = fitz.Rect(0, 330, 720, 648)
    page.show_pdf_page(clip, old_v, 0, clip=clip)
    c.label_pdf(page, old_v[0], clip, clip, 'C  Completeness', 'B')
    c.label_pdf(page, old_v[0], clip, clip, 'D  Mixture', 'C')
    doc.set_metadata({'title': 'Scientific figure', 'author': '[AUTHOR]', 'subject': 'Composition of saved vector panels', 'creator': '[AUTHOR]', 'producer': '[AUTHOR]'})
    path = output / 'figure04-validation-composed.pdf'
    doc.save(path, garbage=4, deflate=True)
    return path


def _without_panel_letters(source, clip, spec):
    """Remove the old letter (under its white patch) and the new letter of one panel title (carried split step)."""
    import fitz
    doc = fitz.open()
    doc.insert_pdf(source, from_page=0, to_page=0)
    page = doc[0]
    found = []
    word = spec['title'].split()[0]
    for block in page.get_text('rawdict')['blocks']:
        for line in block.get('lines', []):
            for span in line['spans']:
                text = ''.join(c['c'] for c in span['chars'])
                if ('  ' + word) in text and clip.contains(fitz.Rect(span['chars'][0]['bbox'])):
                    found.append(fitz.Rect(span['chars'][0]['bbox']))
                elif text == spec['old_panel'] and clip.contains(fitz.Rect(span['bbox'])):
                    found.append(fitz.Rect(span['bbox']))
    assert len(found) == 2, (spec['name'], found)
    page.add_redact_annot(found[0] | found[1], fill=(1, 1, 1))
    page.apply_redactions(images=0, graphics=0, text=0)
    return doc


def split_panel(composed, spec, output):
    """One panel of the composed page, scaled to 720 pt wide, without its letter: PDF, 300 dpi PNG and TIFF."""
    import fitz
    from PIL import Image, PngImagePlugin
    original = fitz.open(composed)
    clip = fitz.Rect(spec['pdf_clip'])
    scale = WIDTH / clip.width
    clean = _without_panel_letters(original, clip, spec)
    doc = fitz.open()
    page = doc.new_page(width=WIDTH, height=clip.height * scale)
    page.show_pdf_page(page.rect, clean, 0, clip=clip)
    doc.set_metadata({'title': spec['name'] + ' - ' + spec['title'], 'author': '[AUTHOR]',
                      'subject': 'Single saved vector panel;scientific content unchanged', 'creator': '[AUTHOR]', 'producer': '[AUTHOR]'})
    path = output / spec['name']
    doc.save(path.with_suffix('.pdf'), garbage=4, deflate=True)
    pix = doc[0].get_pixmap(matrix=fitz.Matrix(DPI / 72, DPI / 72), alpha=False)
    image = Image.frombytes('RGB', (pix.width, pix.height), pix.samples)
    meta = PngImagePlugin.PngInfo()
    meta.add_text('Author', '[AUTHOR]')
    image.save(path.with_suffix('.png'), dpi=(DPI, DPI), pnginfo=meta)
    image.save(path.with_suffix('.tiff'), dpi=(DPI, DPI), compression='tiff_lzw')
    return [path.with_suffix(s).name for s in ['.pdf', '.png', '.tiff']]


def entities(pipeline):
    """Entities wholly absent from the final result (carried from figures/make_figures.py:primary_entities)."""
    plt, one, value, record, save = pipeline.plt, pipeline.one, pipeline.value, pipeline.record, pipeline.save
    ORANGE, BLUE = pipeline.ORANGE, pipeline.BLUE
    rows = pipeline.read('r4a_figure_entity_species.tsv')
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.4, 5.4))
    levels = ['precursor', 'peptide', 'protein_group']
    labels = ['Precursor keys', 'Stripped peptides', 'Protein groups']
    for i, level in enumerate(levels):
        absent = one(rows, filter='F1', level=level, state='wholly_absent')
        count = int(absent['observations'])
        a.barh(i, count, color=ORANGE, height=0.5)
        a.text(count + 40, i, f'{count:,}', va='center')
        record('2', level + ' wholly absent', absent, ['observations'])
        for delta, state, color, label in [(-0.14, 'still_represented', BLUE, 'Still represented'), (0.14, 'wholly_absent', ORANGE, 'Wholly absent')]:
            row = one(rows, filter='F1', level=level, state=state)
            v = value(row, 'observed_yeast_pct')
            b.barh(i + delta, v, color=color, height=0.25, label=label if i == 0 else None)
            b.text(v + 0.008, i + delta, f'{v:.3f}', va='center', fontsize=14)
            record('2', level + ' ' + state, row, ['yeast', 'species_exclusive_denominator', 'observed_yeast_pct', 'absent_vs_still_yeast_fraction_ratio'])
    a.set(yticks=range(3), yticklabels=labels, ylim=(2.6, -0.6), xlabel='Wholly absent cohort entities', title='A  F1 accepted-entity losses')
    a.set_xlim(0, a.get_xlim()[1] * 1.2)
    b.set(yticks=range(3), yticklabels=labels, ylim=(2.6, -0.6), xlabel='Observed yeast fraction (%)', title='B  F1 species-exclusive composition')
    b.set_xlim(0, b.get_xlim()[1] * 1.35)
    b.legend(frameon=False, loc='upper left', bbox_to_anchor=(0, -0.18), fontsize=12)
    save(fig, 'figure02_entities')


def make(single_cell_root, pg250_primary, output):
    """Write FigureS1, FigureS2 and FigureS3 (PDF, PNG, TIFF) and figure-s-data.json into output."""
    import matplotlib
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pg250_primary = Path(pg250_primary)
    data = validation_data(pg250_primary / 'export', pg250_primary / 'fold-change')
    # Each renderer ran in its own process with its own matplotlib settings; rc_context keeps them apart here, and
    # rc_file_defaults gives the settings of a fresh process.
    with matplotlib.rc_context():
        matplotlib.rc_file_defaults()
        render_validation(data, output / 'figure04-validation')
    composed = compose_validation(output)
    written = []
    for spec in SPLITS:
        written += split_panel(composed, spec, output)
    with matplotlib.rc_context():
        matplotlib.rc_file_defaults()
        import importlib
        from mbraudit.figures import panels
        panels = importlib.reload(panels)  # its module-level settings apply inside this context
        staging = output / 'panel-tables'
        staging.mkdir(exist_ok=True)
        shutil.copyfile(Path(single_cell_root) / 'on' / 'accounting' / 'figure_entity_species.tsv', staging / 'r4a_figure_entity_species.tsv')
        panels.TABLES, panels.OUTPUT = staging, output
        entities(panels)
    for suffix in ['.pdf', '.png', '.tiff']:
        shutil.copyfile(output / ('figure02_entities' + suffix), output / ('FigureS2' + suffix))
        written.append('FigureS2' + suffix)
    plotted = {'FigureS1': data['B'], 'FigureS3': data['C'], 'FigureS2': [v for v in panels.VALUES if v['figure'] == '2']}
    (output / 'figure-s-data.json').write_text(json.dumps(plotted, indent=2) + '\n', encoding='utf-8')
    return sorted(written)
