"""Figure 1: panel data from the pipeline outputs, the carried renderers, and the composition of panels A to C.

panel_data() collects every plotted value from the single-cell and 250 pg outputs (no statistic is computed
here). Panels A and B are drawn by the carried headline.render, the completeness panel by the carried
panels.primary_flow, and compose() places them in one figure as the original composition step did.
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
        acc = _tsv(single_cell_root / conf / 'accounting' / 'accounting.tsv').set_index('filter').loc['F3']
        cov = _tsv(single_cell_root / conf / 'accounting' / 'represented_coverage.tsv')
        cov = cov[cov['filter'].eq('F3') & cov.level.eq('precursor')].iloc[0]
        a.append({'configuration': NAMES[conf], 'label': LABELS_A[conf], 'filter': 'F3', 'counts': {
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
        sel = dict(endpoint='T1', mask='F3', scope='pooled', group='all', numerator=numerator, reference='R')
        p = points.loc[(points[list(sel)] == pd.Series(sel)).all(axis=1)].iloc[0]
        lim = limits.loc[(limits[list(sel)] == pd.Series(sel)).all(axis=1) & limits.scheme.eq('two_way') & limits.statistic.eq('ratio')].iloc[0]
        contrasts.append({'contrast': name, 'scheme': 'two_way', 'point': repr(float(p.fraction_ratio)),
                          'low95': repr(float(lim.low95)), 'high95': repr(float(lim.high95))})
    b.append({'configuration': 'PXD056832', 'label': 'PXD056832, timsTOF exact negatives (pooled F3)', 'contrasts': contrasts})
    return {'A': a, 'B': b}


def compose(output):
    """Panels A and B from the headline figure, panel C from the completeness panel (carried composition rule)."""
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
    s = c.node('svg', width=str(hw) + 'pt', height=str(total) + 'pt', viewBox=f'0 0 {hw} {total}', version='1.1')
    s.append(c.node('rect', width=hw, height=total, fill='white'))
    s.append(c.svg_component(h, ['axes_1', 'axes_2', 'axes_3'], 'headline-', 0, 0, hw, hh, h.get('viewBox')))
    s.append(c.svg_component(k, ['axes_4'], 'completeness-', 0, hh + gap, hw, (ch - cut) * scale, f'0 {cut} {cw} {ch - cut}', [('B  Own-row', 'C  Own-row')]))
    old_h = fitz.open(output / 'figure01-headline-v3.pdf')
    old_c = fitz.open(output / 'figure01_accounting.pdf')
    doc = fitz.open()
    page = doc.new_page(width=hw, height=total)
    page.show_pdf_page(fitz.Rect(0, 0, hw, hh), old_h, 0)
    plot_cut = 390.0
    clip = fitz.Rect(0, plot_cut, cw, ch)
    target = fitz.Rect(0, hh + gap + (plot_cut - cut) * scale, hw, total)
    page.show_pdf_page(target, old_c, 0, clip=clip)
    spans = [sp for bl in old_c[0].get_text('dict')['blocks'] if 'lines' in bl for li in bl['lines'] for sp in li['spans'] if sp['text'].startswith('B  Own-row')]
    assert len(spans) == 1
    title = spans[0]
    page.insert_font(fontname='panel-title', fontbuffer=_P(findfont('DejaVu Sans')).read_bytes())
    page.insert_text((title['origin'][0] * scale, hh + gap + (title['origin'][1] - cut) * scale), 'C' + title['text'][1:], fontsize=title['size'] * scale, fontname='panel-title', color=(0, 0, 0))
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
