"""Compose Figure 1 from the saved vector panels (carried from figures/compose_main.py; OUTPUT set by the caller)."""
from pathlib import Path
import copy, hashlib, json, re
from lxml import etree as E
import fitz
from PIL import Image, PngImagePlugin
from matplotlib.font_manager import findfont
SVG = 'http://www.w3.org/2000/svg'
NS = {'s': SVG}
OUTPUT = None


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def node(tag, **kw):
    return E.Element('{' + SVG + '}' + tag, **{k: str(v) for k, v in kw.items()})


def load(n):
    return E.parse(str(OUTPUT / n)).getroot()


def rename_ids(tree, prefix):
    ids = {x.get('id'): prefix + x.get('id') for x in tree.iter() if x.get('id')}
    for x in tree.iter():
        for k, v in list(x.attrib.items()):
            if k == 'id':
                x.set(k, ids[v])
            else:
                for old, new in ids.items():
                    v = v.replace('url(#' + old + ')', 'url(#' + new + ')')
                if v.startswith('#') and v[1:] in ids:
                    v = '#' + ids[v[1:]]
                x.set(k, v)


def svg_component(source, groups, prefix, x, y, width, height, viewbox, relabels=()):
    out = node('svg', x=x, y=y, width=width, height=height, viewBox=viewbox, overflow='hidden')
    for d in source.findall('s:defs', NS):
        out.append(copy.deepcopy(d))
    figure = source.find('s:g', NS)
    for name in groups:
        g = copy.deepcopy(figure.find("s:g[@id='" + name + "']", NS))
        for text in g.findall('.//s:text', NS):
            for old, new in relabels:
                if text.text and text.text.startswith(old):
                    text.text = new + text.text[len(old):]
        out.append(g)
    rename_ids(out, prefix)
    return out


def label_pdf(page, source, source_rect, target_rect, old, new):
    spans = [s for b in source.get_text('rawdict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if ''.join((c['c'] for c in s['chars'])).startswith(old)]
    assert len(spans) == 1, old
    s = spans[0]
    c = s['chars'][0]
    scale = target_rect.width / source_rect.width
    assert abs(scale - target_rect.height / source_rect.height) < 1e-06

    def point(p):
        return fitz.Point(target_rect.x0 + (p[0] - source_rect.x0) * scale, target_rect.y0 + (p[1] - source_rect.y0) * scale)
    a = point(c['bbox'][:2])
    b = point(c['bbox'][2:])
    page.draw_rect(fitz.Rect(a, b), color=None, fill=(1, 1, 1), overlay=True)
    page.insert_font(fontname='panel-label', fontbuffer=Path(findfont('DejaVu Sans')).read_bytes())
    page.insert_text(point(c['origin']), new, fontsize=s['size'] * scale, fontname='panel-label', color=(0, 0, 0), overlay=True)


def save(name, svg, doc):
    p = OUTPUT / name
    p.with_suffix('.svg').write_bytes(E.tostring(svg, xml_declaration=True, encoding='UTF-8') + b'\n')
    doc.set_metadata({'title': 'Scientific figure', 'author': '[AUTHOR]', 'subject': 'Composition of saved vector panels', 'creator': '[AUTHOR]', 'producer': '[AUTHOR]'})
    doc.save(p.with_suffix('.pdf'), garbage=4, deflate=True)
    pix = doc[0].get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False)
    im = Image.frombytes('RGB', [pix.width, pix.height], pix.samples)
    meta = PngImagePlugin.PngInfo()
    meta.add_text('Author', '[AUTHOR]')
    im.save(p.with_suffix('.png'), dpi=(300, 300), pnginfo=meta)
    im.save(p.with_suffix('.tiff'), dpi=(300, 300), compression='tiff_lzw')
    return {p.with_suffix(s).name: sha(p.with_suffix(s)) for s in ['.svg', '.pdf', '.png', '.tiff']}
