"""Species census of a library: label, stripped and predicted are carried from analysis/candidate_ratios.py.

predicted() decodes the binary predicted library PRED into one row per precursor and caches it under ROOT;
both module globals are set by the entry script (plumbing).
"""
import json
import mmap
import struct
from collections import Counter
import numpy as np
import pandas as pd
from mbraudit import species as COMMON
from mbraudit.library import candidate_library as PG
from mbraudit.library import candidate_library as BIN
PRED = None
ROOT = None


def stripped(s):
    out, depth = ([], 0)
    for ch in s:
        if ch in '([':
            depth += 1
        elif ch in ')]':
            depth = max(0, depth - 1)
        elif depth == 0 and 'A' <= ch <= 'Z':
            out.append(ch)
    return ''.join(out)


def predicted():
    cache = ROOT / 'cache/predicted-candidates.parquet'
    if cache.exists():
        return pd.read_parquet(cache)
    maps = PG.load_pg(PRED)
    mapping = maps['idx']
    rows, flags = ([], Counter())
    with PRED.open('rb') as h, mmap.mmap(h.fileno(), 0, access=mmap.ACCESS_READ) as b:
        cursor = BIN.Cur(b)
        header = BIN.parse_header(cursor)
        count = cursor.i32()
        assert count == 1745123 == maps['mapped']
        p = cursor.p
        for i in range(count):
            rid, ch, length = struct.unpack_from('<iii', b, p)
            p += 36 + (4 if header[0] <= -8 else 0)
            flag = struct.unpack_from('<i', b, p)[0]
            p += 8
            n = struct.unpack_from('<i', b, p)[0]
            p += 4
            key = b[p:p + n].decode('utf-8')
            p += n + 16
            nf = struct.unpack_from('<i', b, p)[0]
            p += 4 + nf * 12
            charge = ch & 65535
            assert rid == i and key.endswith(str(charge)) and (2 <= charge <= 4)
            seq = key[:-len(str(charge))]
            aa = stripped(seq)
            assert len(aa) == length & 65535
            rows.append((seq, charge, aa, mapping[i][0], flag))
            flags[flag] += 1
            if i % 100000 == 0:
                print('Parsed', i, flush=True)
        assert len(b) - p == 4 * count + 12, (p, len(b))
        assert struct.unpack_from('<i', b, p)[0] == count
        assert struct.unpack_from('<ii', b, len(b) - 8) == (0, 0)
    df = pd.DataFrame(rows, columns=['Modified.Sequence', 'Precursor.Charge', 'Stripped.Sequence', 'Protein.Ids', 'native_flags'])
    assert not df.duplicated(['Modified.Sequence', 'Precursor.Charge']).any()
    df.to_parquet(cache)
    (ROOT / 'cache/predicted-layout.json').write_text(json.dumps({'records': count, 'mapped': len(mapping), 'payload_end': p, 'file_bytes': PRED.stat().st_size, 'trailer_bytes': 4 * count + 12, 'flags': dict(flags), 'explicit_decoy_names': int(df['Modified.Sequence'].str.contains('DECOY').sum()), 'protein_records': maps['proteinCount'], 'peptide_groups': maps['pgCount']}, indent=2))
    return df


def label(df, registry, role):
    unique = df[['Protein.Ids', 'Stripped.Sequence']].drop_duplicates().copy()
    labels = []
    for start in range(0, len(unique), 1000):
        block = unique.iloc[start:start + 1000]
        labels.extend(registry.labels(block['Protein.Ids'], block['Stripped.Sequence']))
        if start % 10000 == 0:
            print(role, 'species labels', start, '/', len(unique), flush=True)
    unique['species'] = labels
    df = df.merge(unique, on=['Protein.Ids', 'Stripped.Sequence'], validate='many_to_one')
    df['species_before_IL'] = [COMMON.classify(x, registry.target, registry.entrap) for x in df['Protein.Ids']]
    return df
