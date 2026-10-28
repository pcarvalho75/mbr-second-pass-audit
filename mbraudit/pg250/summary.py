"""Whole-pass sums (first pass = retained + lost; final = retained + added) from class counts (carried from analysis/v5_summary.py)."""
import csv
from decimal import Decimal, localcontext
from pathlib import Path


def accepted_passes(rows, organism, mask='F1'):
    """Retained plus lost is first pass;retained plus added is the final pass."""
    selected = {}
    for category in ['R', 'T', 'L']:
        found = [r for r in rows if r['endpoint'] == 'T1' and r['mask'] == mask and r['scope'] == 'organism'
                 and r['group'] == organism and r['tercile'] == 'all' and r.get('context', r.get('class')) == category]
        if len(found) != 1:
            raise ValueError('Expected one organism/class count row: ' + organism + '/' + category)
        selected[category] = found[0]
    result = []
    for phase, classes in [('first', ['R', 'L']), ('final', ['R', 'T'])]:
        absent = sum(int(selected[c]['A']) for c in classes)
        present = sum(int(selected[c]['P']) for c in classes)
        with localcontext() as context:
            context.prec = 50
            fraction = Decimal(absent) / Decimal(absent + present)
        result.append({'organism': organism, 'pass': phase, 'A': absent, 'P': present, 'n': absent + present,
                       'observed_wrong_organism_fraction': str(fraction), 'percent': f'{fraction*100:.3f}'})
    return result


def read(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def added_count_comparison(mixed_rows, own_rows, organism='H', mask='F1'):
    def one(rows):
        found = [r for r in rows if r['endpoint'] == 'T1' and r['mask'] == mask and r['scope'] == 'organism'
                 and r['group'] == organism and r['tercile'] == 'all' and r.get('context', r.get('class')) == 'T']
        if len(found) != 1:
            raise ValueError('Expected one added-class count for ' + organism)
        return found[0]
    mixed, own = one(mixed_rows), one(own_rows)
    a, b = int(mixed['A']), int(own['A'])
    reduction = (Decimal(a - b) / Decimal(a) * 100) if a else None
    return {'organism': organism, 'wrong_added_mixed': a, 'wrong_added_own': b,
            'right_added_mixed': int(mixed['P']), 'right_added_own': int(own['P']),
            'wrong_added_count_reduction_percent': str(reduction) if reduction is not None else 'NA'}
