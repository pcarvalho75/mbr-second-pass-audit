"""Report acceptance, precursor identity and species labels for the single-cell dataset.

Carried from the old analysis/common.py. TARGET_FASTA and ENTRAP_FASTA come from the configuration.
"""
from __future__ import annotations
import hashlib
import os
import numpy as np
import pandas as pd
from scipy.stats import beta,chi2
from mbraudit.config import optional_input
TARGET_FASTA=optional_input('single_cell.target_fasta')
ENTRAP_FASTA=optional_input('single_cell.entrapment_fasta')


ENTRAP_RATIO=0.276123


Q_ACCEPT=0.01


PXD_RUN_PREFIX='20230929_yzc_'


RUN_EXTENSIONS=('.raw.dia','.raw','.mzML','.dia','.d')


REPORT_KEY_COLUMNS=['Run','Modified.Sequence','Precursor.Charge','Q.Value']


def run_name(value: str) -> str:
    base = os.path.basename(str(value).replace('\\', '/'))
    for ext in RUN_EXTENSIONS:
        if base.lower().endswith(ext.lower()):
            return base[:-len(ext)]
    return base


def load_registry(path: str) -> tuple[set[str], dict[str, str]]:
    accs, seqs, acc, buf = (set(), {}, None, [])
    with open(path, encoding='utf8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if acc is not None:
                    seqs[acc] = ''.join(buf)
                head = line[1:].split()[0]
                parts = head.split('|')
                acc = parts[1] if len(parts) >= 2 else head
                accs.add(acc)
                buf = []
            else:
                buf.append(line.strip())
    if acc is not None:
        seqs[acc] = ''.join(buf)
    return (accs, seqs)


def classify(protein_ids: str, target: set[str], entrap: set[str]) -> str:
    ids = [t.strip() for t in str(protein_ids).split(';') if t and t.strip() and (t.strip().lower() != 'nan')]
    if not ids:
        return 'unknown'
    in_t = [i in target for i in ids]
    in_e = [i in entrap and i not in target for i in ids]
    if not all((a or b for a, b in zip(in_t, in_e))):
        return 'unknown'
    if all(in_t):
        return 'target'
    if all(in_e):
        return 'entrapment'
    return 'shared'


class SpeciesRegistry:

    def __init__(self, fasta_target: str=TARGET_FASTA, fasta_entrap: str=ENTRAP_FASTA):
        self.target, seqs = load_registry(fasta_target)
        self.entrap, _ = load_registry(fasta_entrap)
        self.entrap -= self.target
        self.target_il = '#'.join((s.replace('I', 'L') for s in seqs.values()))

    def labels(self, protein_ids, stripped) -> list[str]:
        sp = [classify(p, self.target, self.entrap) for p in protein_ids]
        return ['shared_IL' if s == 'entrapment' and str(st).replace('I', 'L') in self.target_il else s for s, st in zip(sp, stripped)]


def combined_factor(r: float=ENTRAP_RATIO) -> float:
    return 1 + 1 / r


def clopper_pearson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return (float(lo), float(hi))


def poisson_combined_bounds(e: int, n: int, r: float=ENTRAP_RATIO) -> tuple[float, float, float]:
    k = combined_factor(r)
    lo = chi2.ppf(0.025, 2 * e) / 2 if e else 0.0
    hi = chi2.ppf(0.975, 2 * (e + 1)) / 2
    return (e * k / n, lo * k / n, hi * k / n)


def accept_rows(d: pd.DataFrame, q: float=Q_ACCEPT, targets_only: bool=False) -> pd.DataFrame:
    d['Run'] = d['Run'].map(run_name).str.replace(PXD_RUN_PREFIX, '', regex=False)
    keep = d['Q.Value'] <= q
    if targets_only:
        keep = keep & (d['Decoy'] == 0)
    d = d[keep].copy()
    d['key'] = list(zip(d['Modified.Sequence'], d['Precursor.Charge'].astype(int)))
    return d
