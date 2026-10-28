"""The fixed precursor and protein-group q-value filters F1, F2 and F3 (carried from analysis/filters.py)."""
from __future__ import annotations
import pandas as pd


Q_CUT=0.01


FILTERS=['F1','F2','F3']


def filter_mask(df: pd.DataFrame, name: str) -> pd.Series:
    if name not in FILTERS:
        raise KeyError(name)
    mask = df['Q.Value'] <= Q_CUT
    if name in ('F1', 'F2'):
        mask = mask & (df['Global.Q.Value'] <= Q_CUT)
        mask = mask & (df['Lib.Q.Value'] <= Q_CUT)
    if name == 'F2':
        mask = mask & (df['PG.Q.Value'] <= Q_CUT)
        mask = mask & (df['Global.PG.Q.Value'] <= Q_CUT)
        mask = mask & (df['Lib.PG.Q.Value'] <= Q_CUT)
    return mask
