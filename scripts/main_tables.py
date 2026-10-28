"""Tables 1 to 4 of the manuscript, from the pipeline outputs.

Table 1: retained, added and lost observations, represented precursors and completeness (single cells, F1).
Table 2: yeast trap fractions by class, two-way intervals, ratios to retained and estimated upper bounds.
Table 3: wrong-organism identifications in the single-organism runs of the 250 pg dataset (all 54 runs).
Table 4: human-only runs analysed with the other 48 runs and on their own.
Outputs: <output_dir>/tables/table1.tsv ... table4.tsv (full precision) and tables.md (printed rounding).
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
    t3 = tables.table3(root / 'pg250' / 'primary' / 'export')
    pd.DataFrame(t3).to_csv(out / 'table3.tsv', sep='\t', index=False, lineterminator='\n')
    md += ['', '## Table 3', '', '| Intended organism | Class | A | P | n | Observed wrong-organism fraction [two-way 95% interval] |', '| --- | --- | --- | --- | --- | --- |']
    for r in t3:
        md.append(f"| {r['organism']} | {r['class']} | {r['A']} | {r['P']} | {r['n']} | {r['fraction']:.4g} [{r['low95']:.4g}, {r['high95']:.4g}] |")
    t4, comparison = tables.table4(root / 'pg250' / 'primary' / 'export', root / 'pg250' / 'within-H')
    pd.DataFrame(t4).to_csv(out / 'table4.tsv', sep='\t', index=False, lineterminator='\n')
    md += ['', '## Table 4', '', '| Class | With all 54 runs: A/n | % | Human-only runs alone: A/n | % | Human identifications (P): with / alone |', '| --- | --- | --- | --- | --- | --- |']
    for r in t4:
        md.append(f"| {r['class']} | {r['with_A']}/{fmt_int(r['with_n'])} | {r['with_percent']:.3f} | {r['alone_A']}/{fmt_int(r['alone_n'])} | "
                  f"{r['alone_percent']:.3f} | {fmt_int(r['with_P'])} / {fmt_int(r['alone_P'])} |")
    md.append(f"Wrong-organism additions: {comparison['wrong_added_mixed']} with the mixtures, {comparison['wrong_added_own']} alone "
              f"({float(comparison['wrong_added_count_reduction_percent']):.0f}% fewer).")
    (out / 'tables.md').write_text('\n'.join(md) + '\n', encoding='utf-8')
    print('Tables 1 to 4:', out, flush=True)


if __name__ == '__main__':
    main()
