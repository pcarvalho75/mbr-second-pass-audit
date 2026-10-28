"""Compare every recomputed number with the value printed in the manuscript (references/published-values.tsv).

A number matches when it agrees within half a unit of the last printed digit (for example 3.28 matches any
value from 3.275 to 3.285); integers, count pairs and yes/no statements must agree exactly.
Writes <output_dir>/key-numbers.tsv and <output_dir>/verification.tsv; exits with status 1 on any mismatch.
"""
import sys
from decimal import Decimal

from _setup import setup, REPO


def matches(value, printed, rule):
    if rule == 'bool':
        return str(value) == printed
    if rule == 'exact':
        text = printed.replace(',', '')
        if any(c in text for c in '/ '):
            return str(value) == text
        return int(value) == int(text)
    text = printed.replace(',', '')
    digits = len(text.split('.')[1]) if '.' in text else 0
    tolerance = Decimal(5) * Decimal(10) ** (-digits - 1)
    return abs(Decimal(repr(float(value))) - Decimal(text)) <= tolerance + Decimal('1e-12')


def main():
    setup(__doc__)
    import pandas as pd
    from mbraudit import config
    from mbraudit.key_numbers import compute
    root = config.output_root()
    computed = pd.DataFrame(compute(root, REPO / config.value('candidate_counts')))
    computed.to_csv(root / 'key-numbers.tsv', sep='\t', index=False, lineterminator='\n')
    values = dict(zip(computed['id'], computed['value']))
    published = pd.read_csv(REPO / 'references' / 'published-values.tsv', sep='\t', dtype=str, keep_default_na=False)
    rows = []
    for _, p in published.iterrows():
        if p['id'] not in values:
            rows.append({**p, 'recomputed': '', 'status': 'not computed'})
            continue
        v = values[p['id']]
        rows.append({**p, 'recomputed': v, 'status': 'match' if matches(v, p['printed'], p['decimals']) else 'MISMATCH'})
    report = pd.DataFrame(rows)
    report.to_csv(root / 'verification.tsv', sep='\t', index=False, lineterminator='\n')
    summary = report.status.value_counts().to_dict()
    print('Verification:', summary, flush=True)
    for _, r in report[report.status != 'match'].iterrows():
        print('  ', r['status'], r['id'], 'printed', r['printed'], 'recomputed', r['recomputed'])
    sys.exit(0 if summary.get('match', 0) == len(report) else 1)


if __name__ == '__main__':
    main()
