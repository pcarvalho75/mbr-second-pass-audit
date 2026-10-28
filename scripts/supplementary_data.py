"""Supplementary Data 1, the published workbook: a README sheet and tables D1 to D11 (one sheet each, with blocks)
from the pipeline outputs, each value written with the digits of the published tables (mbraudit.workbook.publish).
The v5.1 tables behind the short Supplement tables S3, S5, S6 and S10 are written at full precision to a separate
check workbook.

With --compare REFERENCE.xlsx, every cell of a reference workbook is checked against the recomputed table at full
precision: rows are matched on their label columns, numbers must agree within half a unit of the reference's last
digit (taken from the cell's number format, or from the printed text), and labels, count pairs and words must agree
exactly. Writes <output_dir>/supplementary-data.xlsx, supplement-table-sources.xlsx and, when comparing,
supplementary-data-comparison.tsv.
"""
import re
from decimal import Decimal

from _setup import setup, REPO

NUMBER = re.compile(r'-?\d+(?:\.\d+)?(?:e-?\d+)?')


def norm(text):
    return re.sub(r'\s+', ' ', str(text)).strip().lower()


def digits_of(cell):
    fmt = cell.number_format or 'General'
    if fmt != 'General' and '.' in fmt:
        return len(fmt.split('.')[1].rstrip('%;@'))
    if fmt != 'General':
        return 0
    text = repr(cell.value)
    return len(text.split('.')[1]) if '.' in text and 'e' not in text else 0


def close(mine, printed, places):
    try:
        value = Decimal(repr(float(mine)))
    except (TypeError, ValueError):
        return False
    return abs(value - Decimal(str(printed))) <= Decimal(5) * Decimal(10) ** (-places - 1) + Decimal('1e-12')


def tokens(text):
    parts = NUMBER.split(str(text))
    return [norm(p) for p in parts], NUMBER.findall(str(text))


def same_text(mine, published):
    a_words, a_numbers = tokens(mine)
    b_words, b_numbers = tokens(published)
    if a_words != b_words or len(a_numbers) != len(b_numbers):
        return False
    for x, y in zip(a_numbers, b_numbers):
        places = len(y.split('.')[1]) if '.' in y else 0
        if not close(x, y, places):
            return False
    return True


def compare_cell(mine, cell):
    published = cell.value
    if published is None:
        return True
    if str(mine) == str(published):
        return True
    if isinstance(published, bool):
        return str(mine) == str(published)
    if isinstance(published, (int, float)):
        if isinstance(mine, str):
            if norm(mine) in ('unavailable', 'nan'):
                return False
            return same_text(mine, repr(published))
        return close(mine, published, digits_of(cell))
    if isinstance(mine, (int, float)):
        return same_text(repr(float(mine)) if isinstance(mine, float) else str(mine), published)
    return same_text(mine, published)


def find_header(rows, columns, start=0):
    want = [norm(c) for c in columns]
    for i in range(start, len(rows)):
        cells = [norm(c.value) if c.value is not None else '' for c in rows[i]]
        if cells[:len(want)] == want:
            return i
    return None


def compare(sheets, reference):
    import openpyxl
    wb = openpyxl.load_workbook(reference)
    report = []
    for name, blocks in sheets.items():
        if isinstance(blocks, Exception):
            report.append({'sheet': name, 'block': '', 'row': '', 'column': '', 'published': '', 'recomputed': '', 'status': 'not produced: ' + str(blocks)})
            continue
        if name not in wb.sheetnames:
            continue  # extra sheets (short Supplement tables) have no counterpart in the published workbook
        rows = list(wb[name].iter_rows())
        start = 0
        from mbraudit.workbook import KEYS
        for (title, frame), nkeys in zip(blocks, KEYS[name]):
            columns = list(frame.columns)
            header = find_header(rows, columns, start)
            if header is None:
                report.append({'sheet': name, 'block': title, 'row': '', 'column': '', 'published': '', 'recomputed': '', 'status': 'header not found'})
                continue
            keycols = columns[:nkeys]
            mine = {tuple(norm(v) for v in r[keycols]): r for _, r in frame.iterrows()}
            i = header + 1
            while i < len(rows) and any(c.value is not None for c in rows[i]):
                cells = rows[i][:len(columns)]
                published = {c: cells[j] for j, c in enumerate(columns)}
                key = tuple(norm(published[c].value) for c in keycols)
                found = mine.get(key)
                if found is None:
                    report.append({'sheet': name, 'block': title, 'row': i + 1, 'column': '/'.join(keycols), 'published': ' | '.join(key), 'recomputed': '', 'status': 'row not recomputed'})
                else:
                    for c in columns:
                        if c in keycols:
                            continue
                        ok = compare_cell(found[c], published[c])
                        report.append({'sheet': name, 'block': title, 'row': i + 1, 'column': c, 'published': published[c].value,
                                       'recomputed': found[c], 'status': 'match' if ok else 'MISMATCH'})
                i += 1
            start = i
    return report


def main():
    from pathlib import Path
    args = setup(__doc__, lambda p: p.add_argument('--compare', type=Path, default=None, help='published workbook to check against'))
    import pandas as pd
    from mbraudit import config
    from mbraudit.workbook import build, publish, write_extra
    root = config.output_root()
    sheets = build(root, REPO / 'design' / 'samples54.tsv')
    cells = publish(sheets, root / 'supplementary-data.xlsx')
    write_extra(sheets, root / 'supplement-table-sources.xlsx')
    print('Supplementary Data 1:', root / 'supplementary-data.xlsx', f'({cells} table cells, each displaying its text)', flush=True)
    if args.compare:
        report = pd.DataFrame(compare(sheets, args.compare))
        report.to_csv(root / 'supplementary-data-comparison.tsv', sep='\t', index=False, lineterminator='\n')
        summary = report.groupby('sheet').status.value_counts().unstack(fill_value=0)
        print(summary.to_string(), flush=True)
        print('Total:', report.status.value_counts().to_dict(), flush=True)


if __name__ == '__main__':
    main()
