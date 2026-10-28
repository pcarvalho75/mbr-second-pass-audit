"""Predicted-library candidate census and fixed species map (carried from validation_v5/library_opportunity_audit.py)."""
import argparse
from collections import Counter
import hashlib
import json
import mmap
from pathlib import Path
import re
import struct
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from mbraudit.library.candidate_library import Cur, read_i32, read_i32_vec, read_str
BITS = {'H': 1, 'Y': 2, 'E': 4}
TAX = {'9606': ('H', 'HUMAN'), '559292': ('Y', 'YEAST'), '83333': ('E', 'ECOLI')}
CATEGORIES = ['H', 'Y', 'E', 'shared', 'shared_IL', 'unknown', 'quarantined']

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def write_json(path, value):
    Path(path).write_bytes((json.dumps(value, indent=2, allow_nan=False) + '\n').encode())

def strip_mods(value):
    out = []
    depth = 0
    for ch in value:
        if ch in '([':
            depth += 1
        elif ch in ')]':
            depth -= 1
            if depth < 0:
                raise ValueError('Unbalanced modification')
        elif not depth and 'A' <= ch <= 'Z':
            out.append(ch)
    if depth:
        raise ValueError('Unbalanced modification')
    return ''.join(out)

def fasta_registry(path, registry_path=None):
    proteins = {}
    species = {}
    accession = None
    sequence = []

    def finish():
        if accession is not None:
            if accession in proteins:
                raise ValueError('Duplicate FASTA accession')
            proteins[accession] = ''.join(sequence)
    with Path(path).open(encoding='utf-8-sig') as f:
        for line in f:
            if line.startswith('>'):
                finish()
                parts = line[1:].split()[0].split('|')
                if len(parts) != 3:
                    raise ValueError('Expected reviewed FASTA accession syntax')
                accession = parts[1]
                sequence = []
                tax = re.search('\\bOX=(\\d+)\\b', line)
                if not tax or tax[1] not in TAX:
                    raise ValueError('Unknown taxonomy')
                sp, suffix = TAX[tax[1]]
                if not parts[2].endswith('_' + suffix):
                    raise ValueError('Conflicting taxonomy/suffix')
                species[accession] = sp
            else:
                sequence.append(line.strip())
    finish()
    if registry_path:
        registry = pd.read_csv(registry_path, sep='\t')
        if registry.accession.duplicated().any():
            raise ValueError('Duplicate registry accession')
        if dict(zip(registry.accession, registry.species)) != species:
            raise ValueError('Registry and FASTA disagree')
    return (proteins, species)

def accession_class(ids, registry):
    names = [x.strip() for x in ids.split(';') if x.strip()]
    if not names or any((x not in registry for x in names)):
        return 'unknown'
    found = {registry[x] for x in names}
    return next(iter(found)) if len(found) == 1 else 'shared'

def exact_il_hits(patterns, proteins, registry):
    import ahocorasick
    automaton = ahocorasick.Automaton(ahocorasick.STORE_INTS)
    for value, index in patterns.items():
        automaton.add_word(value, index)
    automaton.make_automaton()
    hits = np.zeros(len(patterns), dtype=np.uint8)
    for accession, sequence in proteins.items():
        bit = BITS[registry[accession]]
        for _, index in automaton.iter(sequence.replace('I', 'L')):
            hits[index] |= bit
    return hits

def primary_label(before, hitmask, conflict=False):
    if conflict:
        return 'quarantined'
    if before in BITS and hitmask & (7 ^ BITS[before]):
        return 'shared_IL'
    return before

def group_map(path, expected, registry):
    mapping = np.full(expected, -1, dtype=np.int32)
    conflict = np.zeros(expected, dtype=bool)
    groups = []
    before = []
    canonical = []
    duplicates = 0
    conflicts = 0
    with Path(path).open('rb') as f:
        sentinel = read_i32(f)
        f.seek(12, 1)
        read_str(f)
        read_str(f)
        nprot = read_i32(f)
        for _ in range(nprot):
            read_i32(f)
            nmods = read_i32(f)
            read_str(f)
            read_str(f)
            read_str(f)
            f.seek(8 + max(0, nmods) * 4, 1)
        npg = read_i32(f)
        for group in range(npg):
            tree = read_i32(f)
            ids = read_str(f)
            read_str(f)
            read_str(f)
            read_i32_vec(f)
            read_i32_vec(f)
            indices = read_i32_vec(f)
            if tree > 0:
                f.seek(tree * 4, 1)
            groups.append(ids)
            before.append(accession_class(ids, registry))
            canonical.append(tuple(sorted(set((x.strip() for x in ids.split(';') if x.strip())))))
            for index in indices:
                if not 0 <= index < expected:
                    raise ValueError('Out-of-range precursor assignment')
                old = int(mapping[index])
                if old >= 0:
                    duplicates += 1
                    if canonical[old] != canonical[group]:
                        if not conflict[index]:
                            conflicts += 1
                        conflict[index] = True
                else:
                    mapping[index] = group
    return (mapping, conflict, groups, before, {'sentinel': sentinel, 'protein_records': nprot, 'peptide_groups': npg, 'mapped': int((mapping >= 0).sum()), 'unmapped': int((mapping < 0).sum()), 'repeat_group_assignments': duplicates, 'conflicting_precursor_assignments': conflicts})

def audit(config):
    path = Path(config['library'])
    out = Path(config['output'])
    out.mkdir(parents=True, exist_ok=True)
    if sha(path) != config['library_sha256']:
        raise ValueError('Library identity mismatch')
    fasta = Path(config['fasta'])
    if sha(fasta) != config['fasta_sha256']:
        raise ValueError('FASTA identity mismatch')
    registry_path = config.get('registry')
    proteins, registry = fasta_registry(fasta, registry_path)
    expected = int(config['expected_precursors'])
    mapping, conflicts, groups, before, layout = group_map(path, expected, registry)
    patterns = {}
    flags = Counter()
    keys = set()
    explicit_decoy = 0
    unmatched = 0
    rawpath = out / 'raw_annotations.parquet'
    writer = None
    rows = []

    def flush():
        nonlocal writer, rows
        table = pa.Table.from_pylist(rows)
        if writer is None:
            writer = pq.ParquetWriter(rawpath, table.schema, compression='zstd')
        writer.write_table(table)
        rows = []
    with path.open('rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as b:
        c = Cur(b)
        version = c.i32()
        c.skip(12)
        libname = c.sstr()
        source_name = c.sstr()
        from mbraudit.library.candidate_library import parse_header
        c = Cur(b)
        header = parse_header(c)
        count = c.i32()
        p = c.p
        if count != expected:
            raise ValueError('Native log and binary count differ')
        for ordinal in range(count):
            rid, ch, length = struct.unpack_from('<iii', b, p)
            p += 36 + (4 if version <= -8 else 0)
            flag = struct.unpack_from('<i', b, p)[0]
            p += 8
            n = struct.unpack_from('<i', b, p)[0]
            p += 4
            key = b[p:p + n].decode('utf-8')
            p += n + 16
            nf = struct.unpack_from('<i', b, p)[0]
            p += 4 + nf * 12
            charge = ch & 65535
            if rid != ordinal or not 2 <= charge <= 4 or (not key.endswith(str(charge))):
                raise ValueError('Unsupported binary record')
            if key in keys:
                raise ValueError('Duplicate exact modified-sequence/charge key')
            keys.add(key)
            seq = key[:-len(str(charge))]
            aa = strip_mods(seq)
            if len(aa) != length & 65535:
                raise ValueError('Length disagrees')
            pattern = aa.replace('I', 'L')
            index = patterns.setdefault(pattern, len(patterns))
            group = int(mapping[ordinal])
            ids = groups[group] if group >= 0 else ''
            prior = before[group] if group >= 0 else 'unknown'
            names = [x.strip() for x in ids.split(';') if x.strip()]
            decoy = any((x.startswith(('DECOY', 'REV_')) for x in names)) or 'DECOY' in seq
            explicit_decoy += int(decoy)
            flags[flag] += 1
            if names and all((x in proteins for x in names)):
                unmatched += int(not any((pattern in proteins[x].replace('I', 'L') for x in names)))
            rows.append({'Modified.Sequence': seq, 'Precursor.Charge': charge, 'Stripped.Sequence': aa, 'Protein.Ids': ids, 'species_before_IL': prior, 'pattern_index': index, 'annotation_conflict': bool(conflicts[ordinal]), 'Decoy': int(decoy), 'native_flags': flag})
            if len(rows) == 50000:
                flush()
            if ordinal % 1000000 == 0:
                print(config['label'], 'decoded', ordinal, flush=True)
        if rows:
            flush()
        writer.close()
        if len(b) - p != 4 * count + 12 or struct.unpack_from('<i', b, p)[0] != count or struct.unpack_from('<ii', b, len(b) - 8) != (0, 0):
            raise ValueError('Unsupported binary trailer')
        layout.update({'version': version, 'reserved_and_flags': list(header[1:4]), 'records': count, 'unique_keys': len(keys), 'payload_end': p, 'trailer_bytes': len(b) - p, 'flags': dict(flags), 'explicit_decoys': explicit_decoy, 'not_in_assigned_proteins_IL': unmatched})
    del keys, mapping, conflicts
    print(config['label'], 'exact I/L matching', len(patterns), 'patterns', flush=True)
    hits = exact_il_hits(patterns, proteins, registry)
    writer = None
    counts = Counter()
    conversions = Counter()
    zero_hits = 0
    labeled = out / 'predicted_species.parquet'
    for batch in pq.ParquetFile(rawpath).iter_batches(batch_size=100000):
        frame = batch.to_pandas()
        labels = []
        for prior, index, conflict, charge, decoy in frame[['species_before_IL', 'pattern_index', 'annotation_conflict', 'Precursor.Charge', 'Decoy']].itertuples(index=False, name=None):
            hit = int(hits[index])
            label = primary_label(prior, hit, conflict)
            labels.append(label)
            zero_hits += int(hit == 0)
            if not decoy:
                for scope in ['all', str(charge)]:
                    counts['accession_only', scope, 'quarantined' if conflict else prior] += 1
                    counts['symmetric_IL_exclusive', scope, label] += 1
                    legacy = 'shared_IL' if prior == 'Y' and hit & 1 else prior
                    counts['historical_yeast_IL', scope, 'quarantined' if conflict else legacy] += 1
                if label == 'shared_IL':
                    conversions[prior] += 1
        frame['species'] = labels
        table = pa.Table.from_pandas(frame.drop(columns='pattern_index'), preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(labeled, table.schema, compression='zstd')
        writer.write_table(table)
    writer.close()
    if unmatched or zero_hits:
        raise ValueError('Peptide/FASTA sequence validation failed')
    result = []
    for policy in ['accession_only', 'symmetric_IL_exclusive', 'historical_yeast_IL']:
        for charge in ['all', '2', '3', '4']:
            row = {'library': config['label'], 'policy': policy, 'charge': charge}
            row.update({name: int(counts[policy, charge, name]) for name in CATEGORIES})
            row['total_target_keys'] = sum((row[name] for name in CATEGORIES))
            h, y, e = (row['H'], row['Y'], row['E'])
            row['r_Y_over_H'] = y / h if h else np.nan
            row['r_E_over_H'] = e / h if h else np.nan
            for sp, present, absent in [('H', h, y + e), ('Y', y, h + e), ('E', e, h + y)]:
                row['opportunity_absent_over_present_for_' + sp] = absent / present if present else np.nan
            result.append(row)
    pd.DataFrame(result).to_csv(out / 'opportunities.tsv', sep='\t', index=False, lineterminator='\n', na_rep='NA')
    receipt = {'analysis': 'predicted-library opportunity audit; no search endpoint', 'library_label': config['label'], 'library_sha256': sha(path), 'fasta_sha256': sha(fasta), 'registry_sha256': sha(registry_path) if registry_path else None, 'producer_sha256': sha(__file__), 'decoder_sha256': sha(Path(__file__).with_name('candidate_library.py')), 'proteins_by_species': dict(Counter(registry.values())), 'layout': layout, 'IL_shared_from_accession_exclusive': dict(conversions), 'distinct_IL_patterns': len(patterns), 'primary_policy': 'symmetric_IL_exclusive; exact precursor keys unchanged; match within proteins', 'historical_policy': 'yeast-only I/L exclusion retained separately for study continuity', 'endpoint_policy': 'opportunity ratios only; exact-negative endpoints never use r', 'annotation_sha256': sha(labeled), 'opportunities_sha256': sha(out / 'opportunities.tsv')}
    write_json(out / 'audit_receipt.json', receipt)
    print(json.dumps(receipt), flush=True)
    return receipt
