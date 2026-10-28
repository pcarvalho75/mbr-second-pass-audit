"""Figure 1: panel data from the pipeline outputs, the carried renderers, and the composition of panels A to C.

panel_data() collects every plotted value from the single-cell and 250 pg outputs (no statistic is computed
here); its keys 'A' and 'B' are the two panels of the headline renderer, the accounting bars (Figure 1A) and the
class-relative fractions (Figure 1C). Panel A (accounting) and panel C (class-relative fractions) are drawn by the
carried headline.render, panel B (completeness) by the carried panels.primary_flow. compose() stacks A, B and C
from top to bottom; it uses the original composition rule (vector crops of the saved panels, same page size) with
the completeness panel moved between the other two, so that the panels follow the order of the text.
"""
import json
import shutil
from pathlib import Path

import pandas as pd

LABELS_A = {'on': 'Original (refinement on)', 'off': '--no-refine-q (refinement off)', 'standard': 'Standard library (refinement off)'}
NAMES = {'on': 'original', 'off': 'no_refine', 'standard': 'standard_library'}


def _tsv(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', dtype={'filter': str})


def panel_data(single_cell_root, pg250_export):
    single_cell_root, pg250_export = Path(single_cell_root), Path(pg250_export)
    a, b = [], []
    for conf in ['on', 'off', 'standard']:
        acc = _tsv(single_cell_root / conf / 'accounting' / 'accounting.tsv').set_index('filter').loc['F1']
        cov = _tsv(single_cell_root / conf / 'accounting' / 'represented_coverage.tsv')
        cov = cov[cov['filter'].eq('F1') & cov.level.eq('precursor')].iloc[0]
        a.append({'configuration': NAMES[conf], 'label': LABELS_A[conf], 'filter': 'F1', 'counts': {
            'first': int(acc.first_observations), 'final': int(acc.final_observations),
            'retained': int(acc.retained_observations), 'second_pass_only': int(acc.second_pass_only_observations),
            'first_pass_only_in_library': int(acc.lost_library_present),
            'first_pass_only_outside_library': int(acc.lost_library_excluded),
            'first_represented_precursors': int(cov.first_represented_entities),
            'final_represented_precursors': int(cov.final_represented_entities)}})
        ratios = _tsv(single_cell_root / conf / 'intervals' / 'class_ratio_intervals.tsv')
        contrasts = []
        for name in ['second_pass_only_vs_retained', 'lost_vs_retained']:
            row = ratios[ratios.contrast.eq(name) & ratios.scheme.eq('two_way')].iloc[0]
            contrasts.append({'contrast': name, 'scheme': 'two_way', 'point': repr(float(row.point)),
                              'low95': repr(float(row.low95)), 'high95': repr(float(row.high95))})
        b.append({'configuration': NAMES[conf], 'label': 'Astral, yeast entrapment: ' + LABELS_A[conf], 'contrasts': contrasts})
    points = pd.read_csv(pg250_export / 'contrasts.tsv', sep='\t', float_precision='round_trip')
    limits = pd.read_csv(pg250_export / 'contrast_intervals.tsv', sep='\t', float_precision='round_trip')
    contrasts = []
    for numerator, name in [('T', 'second_pass_only_vs_retained'), ('L', 'lost_vs_retained')]:
        sel = dict(endpoint='T1', mask='F1', scope='pooled', group='all', numerator=numerator, reference='R')
        p = points.loc[(points[list(sel)] == pd.Series(sel)).all(axis=1)].iloc[0]
        lim = limits.loc[(limits[list(sel)] == pd.Series(sel)).all(axis=1) & limits.scheme.eq('two_way') & limits.statistic.eq('ratio')].iloc[0]
        contrasts.append({'contrast': name, 'scheme': 'two_way', 'point': repr(float(p.fraction_ratio)),
                          'low95': repr(float(lim.low95)), 'high95': repr(float(lim.high95))})
    b.append({'configuration': 'PXD056832', 'label': 'PXD056832, timsTOF exact negatives (pooled F1)', 'contrasts': contrasts})
    return {'A': a, 'B': b}


def compose(output):
    """Panel A from the top of the headline figure, panel B from the completeness panel, panel C from the bottom of
    the headline figure (carried composition rule; only the vertical order changes, the page size is unchanged)."""
    import fitz
    from pathlib import Path as _P
    from matplotlib.font_manager import findfont
    from mbraudit.figures import compose as c
    c.OUTPUT = output
    h = c.load('figure01-headline-v3.svg')
    k = c.load('figure01_accounting.svg')
    hw, hh = map(float, h.get('viewBox').split()[2:])
    cw, ch = map(float, k.get('viewBox').split()[2:])
    cut = 360.0
    gap = 16.0
    scale = hw / cw
    total = hh + gap + (ch - cut) * scale
    old_h = fitz.open(output / 'figure01-headline-v3.pdf')
    old_c = fitz.open(output / 'figure01_accounting.pdf')
    # The headline figure keeps the placement of the original composition (its whole page fitted into
    # Rect(0, 0, hw, hh), as page.show_pdf_page does) and is split in the white band between the configuration
    # labels under panel A and the title of the class-relative panel.
    source = old_h[0].rect
    f = min(hw / source.width, hh / source.height)
    cx, cy = hw / 2 - source.width / 2 * f, hh / 2 - source.height / 2 * f
    spans = [sp for bl in old_h[0].get_text('dict')['blocks'] if 'lines' in bl for li in bl['lines'] for sp in li['spans']]
    above = max(sp['bbox'][3] for sp in spans if sp['text'] in ('(refinement on)', '(refinement off)'))
    below = [sp['bbox'][1] for sp in spans if sp['text'].startswith('C  Class-relative')]
    assert len(below) == 1 and below[0] > above
    cut_h = (above + below[0]) / 2
    split = cy + cut_h * f
    shift = gap + (ch - cut) * scale
    panel_b = split + gap
    panel_c = split + shift
    s = c.node('svg', width=str(hw) + 'pt', height=str(total) + 'pt', viewBox=f'0 0 {hw} {total}', version='1.1')
    s.append(c.node('rect', width=hw, height=total, fill='white'))
    s.append(c.svg_component(h, ['axes_1'], 'headline-', 0, 0, hw, split, f'0 0 {hw} {split}'))
    s.append(c.svg_component(k, ['axes_4'], 'completeness-', 0, panel_b, hw, (ch - cut) * scale, f'0 {cut} {cw} {ch - cut}'))
    s.append(c.svg_component(h, ['axes_2', 'axes_3'], 'headline-lower-', 0, panel_c, hw, hh - split, f'0 {split} {hw} {hh - split}'))
    doc = fitz.open()
    page = doc.new_page(width=hw, height=total)
    page.show_pdf_page(fitz.Rect(cx, cy, cx + source.width * f, split), old_h, 0, clip=fitz.Rect(0, 0, source.width, cut_h))
    plot_cut = 390.0
    clip = fitz.Rect(0, plot_cut, cw, ch)
    target = fitz.Rect(0, panel_b + (plot_cut - cut) * scale, hw, panel_c)
    page.show_pdf_page(target, old_c, 0, clip=clip)
    spans = [sp for bl in old_c[0].get_text('dict')['blocks'] if 'lines' in bl for li in bl['lines'] for sp in li['spans'] if sp['text'].startswith('B  Own-row')]
    assert len(spans) == 1
    title = spans[0]
    page.insert_font(fontname='panel-title', fontbuffer=_P(findfont('DejaVu Sans')).read_bytes())
    page.insert_text((title['origin'][0] * scale, panel_b + (title['origin'][1] - cut) * scale), title['text'], fontsize=title['size'] * scale, fontname='panel-title', color=(0, 0, 0))
    page.show_pdf_page(fitz.Rect(cx, panel_c, cx + source.width * f, cy + source.height * f + shift), old_h, 0, clip=fitz.Rect(0, cut_h, source.width, source.height))
    return c.save('figure1', s, doc)


def make(single_cell_root, pg250_export, output):
    from mbraudit.figures import headline, panels
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    data = panel_data(single_cell_root, pg250_export)
    (output / 'figure1-data.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    staging = output / 'panel-tables'
    staging.mkdir(exist_ok=True)
    source = Path(single_cell_root) / 'on' / 'accounting'
    for name in ['figure_flow.tsv', 'central_accounting.tsv', 'figure_completeness.tsv']:
        shutil.copyfile(source / name, staging / ('r4a_' + name))
    panels.TABLES, panels.OUTPUT = staging, output
    panels.save(headline.render(panels, data), 'figure01-headline-v3')
    panels.primary_flow()
    return compose(output)
