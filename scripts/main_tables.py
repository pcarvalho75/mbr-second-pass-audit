"""Tables 1 and 2 of the manuscript, from the pipeline outputs.

Table 1: retained, added and lost observations, represented precursors and completeness (single cells, F1).
Table 2: yeast trap fractions by class, two-way intervals, ratios to retained and estimated upper bounds.
The 250 pg rows that were Tables 3 and 4 of the earlier manuscript are now Figure 2 (scripts/figure2.py, every value
in figure2-data.tsv) and Table D12 of Supplementary Data 1 (scripts/supplementary_data.py).
Outputs: <output_dir>/tables/table1.tsv, table2.tsv, table2-final-result.tsv (full precision) and tables.md
(printed rounding).
"""
from _setup import setup, REPO


def fmt_int(v):
    return f'{int(v):,}'


def main():
    setup(__doc__)
    import pandas as pd
    from mbraudit import config, tables
    root = config.output_root()
    out = config.output_dir('tables')
    md = []
    t1 = tables.table1(root / 'single-cell')
    pd.DataFrame(t1).to_csv(out / 'table1.tsv', sep='\t', index=False, lineterminator='\n')
    md += ['## Table 1', '', '| Measure | Refinement on | Refinement off | Standard library, refinement off |', '| --- | --- | --- | --- |']
    for r in t1:
        cells = [f'{r[c]:.2f}' if '%' in r['measure'] else fmt_int(r[c]) for c in tables.CONFIGS]
        md.append('| ' + ' | '.join([r['measure']] + cells) + ' |')
    t2, whole = tables.table2(root / 'single-cell', REPO / config.value('candidate_counts'))
    pd.DataFrame(t2).to_csv(out / 'table2.tsv', sep='\t', index=False, lineterminator='\n')
    pd.DataFrame(whole).to_csv(out / 'table2-final-result.tsv', sep='\t', index=False, lineterminator='\n')
    md += ['', '## Table 2', '', '| Class / setting | Yeast/n | Observed (lower bound), % | Two-way 95% interval, % | Ratio to retained | Combined estimate (upper bound), % |',
           '| --- | --- | --- | --- | --- | --- |']
    for r in t2:
        md.append(f"| {r['class']} / {r['configuration']} | {r['yeast']}/{r['n']} | {r['observed_percent']:.3f} | "
                  f"{r['interval_low_percent']:.3f} to {r['interval_high_percent']:.3f} | {r['ratio_to_retained']:.3f} | {r['upper_bound_percent']:.3f} |")
    for w in whole:
        md.append(f"Final result as a whole, {w['configuration']}: {w['yeast']}/{w['denominator']}, observed {w['lower_percent']}%, upper bound {w['estimated_upper_percent']}% (r = {w['r']:.4f})")
    (out / 'tables.md').write_text('\n'.join(md) + '\n', encoding='utf-8')
    print('Tables 1 and 2:', out, flush=True)


if __name__ == '__main__':
    main()
