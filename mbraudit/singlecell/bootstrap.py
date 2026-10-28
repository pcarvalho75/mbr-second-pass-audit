"""Recipient (cell) and precursor-block resampling with exact multiplicities (carried from analysis/block_bootstrap.py)."""
import numpy as np


SEED = 20261018


DRAWS = 2000


FACTOR = 1 + 1 / 0.276123


def multiplicities(rng, units, draws):
    indices = rng.integers(0, units, size=(draws, units))
    return np.stack([np.bincount(row, minlength=units) for row in indices]).astype(float)


def ratio(m_e, m_n, h_e, h_n):
    m_e, m_n, h_e, h_n = np.broadcast_arrays(m_e, m_n, h_e, h_n)
    result = np.full(m_e.shape, np.nan, dtype=float)
    valid = (m_n > 0) & (h_n > 0) & (h_e > 0)
    np.divide(m_e * h_n, m_n * h_e, out=result, where=valid)
    result[(m_n > 0) & (h_n > 0) & (h_e == 0) & (m_e > 0)] = np.inf
    return result


def fdp(e, n):
    e, n = np.broadcast_arrays(e, n)
    result = np.full(e.shape, np.nan, dtype=float)
    np.divide(e * FACTOR, n, out=result, where=n > 0)
    return result


def interval(values):
    values = np.asarray(values, dtype=float)
    undefined = int(np.isnan(values).sum())
    infinite = int(np.isposinf(values).sum())
    result = {'draws': len(values), 'undefined_draws': undefined, 'infinite_draws': infinite, 'low95': np.nan, 'high95': np.nan, 'interval_status': 'not estimable: undefined draw' if undefined else 'estimated'}
    if undefined:
        return result
    ordered = np.sort(values)

    def quantile(p):
        index = (len(ordered) - 1) * p
        lo = int(np.floor(index))
        hi = int(np.ceil(index))
        weight = index - lo
        if lo == hi or weight == 0:
            return float(ordered[lo])
        if ordered[lo] == ordered[hi]:
            return float(ordered[lo])
        if np.isposinf(ordered[hi]):
            return np.inf
        return float((1 - weight) * ordered[lo] + weight * ordered[hi])
    result.update(low95=quantile(0.025), high95=quantile(0.975))
    return result


def matrices(frames, cells):
    keys = sorted(set().union(*(set(frame.key) for frame in frames)))
    ci = {cell: i for i, cell in enumerate(cells)}
    ki = {key: i for i, key in enumerate(keys)}
    result = []
    for frame in frames:
        indexes = ([ci[c] for c in frame.Run], [ki[k] for k in frame.key])
        n = np.zeros((len(cells), len(keys)), dtype=float)
        e = n.copy()
        d = n.copy()
        np.add.at(n, indexes, 1)
        np.add.at(e, indexes, frame.species.eq('entrapment').to_numpy(dtype=float))
        np.add.at(d, indexes, frame.species.isin(['target', 'entrapment']).to_numpy(dtype=float))
        result.append((n, e, d))
    return (result, keys)


def bootstrap(frames, cells, kind, scheme, *, draws=DRAWS, seed=SEED, batch=100, checkpoint=None):
    arrays, keys = matrices(frames, cells)
    if not keys:
        return np.full(draws, np.nan)
    rng = np.random.default_rng(seed)
    cw = multiplicities(rng, len(cells), draws) if scheme in ['recipient', 'two_way'] else None
    result = []
    for start in range(0, draws, batch):
        stop = min(draws, start + batch)
        pw = multiplicities(rng, len(keys), stop - start) if scheme in ['precursor', 'two_way'] else None
        totals = []
        for channels in arrays:
            part = []
            for channel in channels:
                if scheme == 'recipient':
                    value = cw[start:stop] @ channel.sum(axis=1)
                elif scheme == 'precursor':
                    value = pw @ channel.sum(axis=0)
                elif scheme == 'two_way':
                    value = (cw[start:stop] @ channel * pw).sum(axis=1)
                else:
                    raise ValueError('Unknown resampling scheme')
                part.append(value)
            totals.append(part)
        if kind == 'composition':
            h, m = totals
            value = ratio(m[1], m[0], h[1], h[0])
        elif kind == 'fdp':
            value = fdp(totals[0][1], totals[0][2])
        else:
            raise ValueError('Unknown statistic')
        result.extend(value.tolist())
        if checkpoint is not None:
            checkpoint(stop, draws)
    return np.asarray(result)
