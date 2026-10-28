"""The fixed precursor and protein-group q-value filters F0, F3 and F4 (carried from analysis/filters.py)."""
from __future__ import annotations
import pandas as pd


Q_CUT=0.01


FILTERS=['F0','F1','F2','F3','F4']


def filter_mask(df: pd.DataFrame, name: str) -> pd.Series:
    mask = df['Q.Value'] <= Q_CUT
    if name in ('F1', 'F3', 'F4'):
        mask = mask & (df['Global.Q.Value'] <= Q_CUT)
    if name in ('F2', 'F3', 'F4'):
        mask = mask & (df['Lib.Q.Value'] <= Q_CUT)
    if name == 'F4':
        mask = mask & (df['PG.Q.Value'] <= Q_CUT)
        mask = mask & (df['Global.PG.Q.Value'] <= Q_CUT)
        mask = mask & (df['Lib.PG.Q.Value'] <= Q_CUT)
    if name not in FILTERS:
        raise KeyError(name)
    return mask
