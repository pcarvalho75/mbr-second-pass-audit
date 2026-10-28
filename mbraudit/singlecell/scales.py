"""Candidate-ratio scaling of observed yeast fractions into estimated upper bounds.

scale_point is carried from analysis/search_space_scales.py and final_single_cell_bound from analysis/v5_summary.py.
"""
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from fractions import Fraction


def scale_point(f, human, yeast, configuration, category):
    """Retain the released 50-digit convention and the exact-rational correction."""
    with localcontext() as context:
        context.prec = 50
        context.rounding = ROUND_HALF_EVEN
        value = Decimal(f) * (1 + Decimal(human) / Decimal(yeast))
        percent = value * 100
        if configuration == 'off' and category == 'lost':
            exact = Fraction(f) * (1 + Fraction(human, yeast))
            value = Decimal(exact.numerator) / Decimal(exact.denominator)
            exact *= 100
            percent = Decimal(exact.numerator) / Decimal(exact.denominator)
        return str(value), str(percent), f'{percent:.3f}'


def final_single_cell_bound(rows, human_candidates, yeast_candidates, mask='F1'):
    selected = [r for r in rows if r['filter'] == mask and r['class'] in ['retained', 'second_pass_only']]
    if len(selected) != 2:
        raise ValueError('Expected retained and added class counts')
    yeast = sum(int(r['yeast']) for r in selected)
    denominator = sum(int(r['species_exclusive_denominator']) for r in selected)
    with localcontext() as context:
        context.prec = 50
        lower = Decimal(yeast) / Decimal(denominator)
        upper = lower * (1 + Decimal(human_candidates) / Decimal(yeast_candidates))
        return {'yeast': yeast, 'denominator': denominator, 'lower_fraction': str(lower),
                'estimated_upper_fraction': str(upper), 'lower_percent': f'{lower*100:.3f}',
                'estimated_upper_percent': f'{upper*100:.3f}'}
