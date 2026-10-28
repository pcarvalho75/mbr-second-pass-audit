"""A small synthetic dataset in the layout of the real inputs, for the end-to-end test.

Writes DIA-NN-like parquet reports (first pass, final result, empirical library) for the 64-cell search and the 16
pool searches, two FASTA files, and the 250 pg inputs (54-run and replication co-searches, three within-organism
searches, a species map, designs and export-status tables), plus a configuration file that points at all of it.
Every name starts with SYN_ where the code allows it; nothing here resembles the real data beyond its layout.

    python tests/synthetic.py OUTPUT_FOLDER
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
QF = ['Q.Value', 'Global.Q.Value', 'Lib.Q.Value', 'PG.Q.Value', 'Global.PG.Q.Value', 'Lib.PG.Q.Value']
AMINO = 'ACDEFGHKLMNPQRSTVWY'


def peptide(rng, length=9):
    return ''.join(rng.choice(list('ACDEFGHMNPQRSTVWY'), length - 1)) + 'K'


def proteins(rng, prefix, n, per=3):
    entries = {}
    for i in range(n):
        peps = [peptide(rng) for _ in range(per)]
        entries[f'{prefix}{i:04d}'] = peps
    return entries


def write_fasta(path, entries, suffix):
    with open(path, 'w', encoding='utf-8', newline='\n') as handle:
        for accession, peps in entries.items():
            handle.write(f'>sp|{accession}|{accession}_{suffix}\n' + 'M' + 'G'.join(peps) + '\n')


def frame(rows):
    return pd.DataFrame(rows, columns=['Run', 'Modified.Sequence', 'Stripped.Sequence', 'Precursor.Charge', 'Decoy', 'Protein.Ids',
                                       'Protein.Group', *QF, 'Precursor.Quantity', 'Precursor.Normalised'])


def row(run, key, accession, q=0.001, lib=0.0, quantity=1000.0, decoy=0):
    return [run, key[0], key[0], key[1], decoy, accession, accession, q, q, lib, q, q, lib, quantity, quantity]


def search(folder, runs, catalogue, rng, presence=0.55, gain=0.5):
    """catalogue: list of (key, accession, weight). Writes first pass, final result and empirical library."""
    folder.mkdir(parents=True, exist_ok=True)
    first, final = [], []
    for run in runs:
        for key, accession, weight in catalogue:
            in_first = rng.random() < presence * weight
            q = 0.02 if rng.random() < 0.05 else 0.001
            quantity = float(rng.lognormal(7, 1)) if rng.random() > 0.03 else 0.0
            if in_first:
                first.append(row(run, key, accession, q=q, lib=0.0, quantity=quantity))
            if (in_first and rng.random() < 0.85) or (not in_first and rng.random() < gain * weight):
                final.append(row(run, key, accession, q=0.001 if rng.random() > 0.05 else 0.02, lib=0.004, quantity=quantity))
        first.append(row(run, ('SYN_DECOY_K', 2), 'SYN_DECOY', decoy=1))
    library = [{'Modified.Sequence': k[0], 'Precursor.Charge': k[1], 'Decoy': 0, 'Q.Value': 0.001, 'Stripped.Sequence': k[0], 'Protein.Ids': a}
               for i, (k, a, _) in enumerate(catalogue) if i % 7]
    frame(first).to_parquet(folder / 'report-first-pass.parquet', index=False)
    frame(final).to_parquet(folder / 'report.parquet', index=False)
    pd.DataFrame(library).to_parquet(folder / 'report-lib.parquet', index=False)


def single_cell(root, rng):
    from mbraudit.singlecell.design import LINES, M6, H16, M16, H6_RUNS
    human, yeast = proteins(rng, 'SYNH', 40), proteins(rng, 'SYNY', 12)
    write_fasta(root / 'human.fasta', human, 'HUMAN')
    write_fasta(root / 'human-yeast.fasta', {**human, **yeast}, 'HUMAN')
    with open(root / 'human-yeast.fasta', 'a', encoding='utf-8', newline='\n') as handle:
        pass
    catalogue = [((p, 2), a, 1.0) for a, peps in human.items() for p in peps] + [((p, 2), a, 0.15) for a, peps in yeast.items() for p in peps]
    cells = [c for members in LINES.values() for c in members]
    for name in ['on', 'off', 'standard']:
        search(root / 'single-cell' / name, cells, catalogue, rng)
    lines = list(LINES.values())
    pools = {alias: lines[i][:6] for i, alias in enumerate(H6_RUNS)}
    pools.update(M6)
    pools.update(H16)
    pools.update(M16)
    for alias, members in pools.items():
        search(root / 'pools' / alias, members, catalogue, rng)
    return list(pools)


DESIGN_COLUMNS = ['sample', 'kind', 'group', 'filename_date', 'present_species']


def design54(date):
    rows = []
    for group in ['5H1Y4E', '5H2Y3E', '4H2Y2E', '5H3Y2E', '5H4Y1E']:
        rows += [{'sample': f'SYN_{date}_{group}_{i}', 'kind': 'mixture', 'group': group, 'filename_date': date, 'present_species': 'H;Y;E'} for i in range(6)]
    for sp in 'HYE':
        rows += [{'sample': f'SYN_single_{sp}_{i}', 'kind': 'single', 'group': sp, 'filename_date': '20250312', 'present_species': sp} for i in range(6)]
    rows += [{'sample': f'SYN_blank_{i}', 'kind': 'blank', 'group': 'blank', 'filename_date': '20250312', 'present_species': ''} for i in range(6)]
    return pd.DataFrame(rows, columns=DESIGN_COLUMNS)


def pg_search(folder, design, keys, rng):
    folder.mkdir(parents=True, exist_ok=True)
    reports = {'first': [], 'final': []}
    for _, d in design.iterrows():
        for (key, sp) in keys:
            if d.kind == 'blank':
                weight = 0.02
            elif d.kind == 'single':
                weight = 0.9 if sp == d.present_species else 0.05
            else:
                weight = 0.9
            share = {'H': 5, 'Y': int(d.group[2]) if d.kind == 'mixture' else 1, 'E': int(d.group[4]) if d.kind == 'mixture' else 1}.get(sp, 5)
            quantity = float(share * rng.lognormal(6, 0.3))
            in_first = rng.random() < weight
            for stage in ['first', 'final']:
                present = in_first if stage == 'first' else (in_first or rng.random() < weight * 0.5)
                if present:
                    q = 0.001 if rng.random() > 0.05 else 0.02
                    reports[stage].append({'Run': d['sample'] + '.d', 'Modified.Sequence': key[0], 'Precursor.Charge': key[1], 'Decoy': 0,
                                           **{f: q for f in QF}, 'Precursor.Quantity': quantity, 'Precursor.Normalised': quantity})
        reports['first'].append({'Run': d['sample'] + '.d', 'Modified.Sequence': 'SYN_DECOY', 'Precursor.Charge': 2, 'Decoy': 1,
                                 **{f: 0.001 for f in QF}, 'Precursor.Quantity': 1.0, 'Precursor.Normalised': 1.0})
    pd.DataFrame(reports['first']).to_parquet(folder / 'report-first-pass.parquet', index=False)
    pd.DataFrame(reports['final']).to_parquet(folder / 'report.parquet', index=False)
    pd.DataFrame([{'Modified.Sequence': k[0], 'Precursor.Charge': k[1], 'Decoy': 0} for k, _ in keys[::2]]).to_parquet(folder / 'report-lib.parquet', index=False)


def pg250(root, rng):
    keys = [((f'SYN{sp}PEPTIDE{i}K', 2), sp) for sp in 'HYE' for i in range(25)]
    keys += [(('SYNSHAREDK', 2), 'shared')]
    species = pd.DataFrame([{'Modified.Sequence': k[0], 'Precursor.Charge': k[1], 'species': sp, 'Stripped.Sequence': k[0],
                             'Protein.Ids': 'SYN_' + sp} for k, sp in keys])
    species.to_parquet(root / 'predicted_species.parquet', index=False)
    for study, date in [('primary', '20250116'), ('replication', '20250109')]:
        design = design54(date)
        design.to_csv(root / f'design-{study}.tsv', sep='\t', index=False, lineterminator='\n')
        pd.DataFrame([{'sample': s, 'stage': st, 'status': 'complete'} for s in design['sample'] for st in ['first', 'final']]).to_csv(
            root / f'exports-{study}.tsv', sep='\t', index=False, lineterminator='\n')
        pg_search(root / study, design, keys, rng)
    for sp in 'HYE':
        design = design54('20250116')
        design = design[(design.kind == 'single') & (design.present_species == sp)]
        design.to_csv(root / f'design-within-{sp}.tsv', sep='\t', index=False, lineterminator='\n')
        pd.DataFrame([{'sample': s, 'stage': st, 'status': 'complete'} for s in design['sample'] for st in ['first', 'final']]).to_csv(
            root / f'exports-within-{sp}.tsv', sep='\t', index=False, lineterminator='\n')
        pg_search(root / ('within-' + sp), design, keys, rng)


def build(root):
    root = Path(root)
    rng = np.random.default_rng(20261028)
    (root / 'pg250').mkdir(parents=True, exist_ok=True)
    pools = single_cell(root, rng)
    pg250(root / 'pg250', rng)
    config = {'output_dir': str(root / 'outputs'),
              'single_cell': {'target_fasta': str(root / 'human.fasta'), 'entrapment_fasta': str(root / 'human-yeast.fasta'),
                              'configurations': {c: str(root / 'single-cell' / c) for c in ['on', 'off', 'standard']},
                              'pools': {alias: str(root / 'pools' / alias) for alias in pools}},
              'candidate_counts': 'references/candidate-counts.tsv',
              'pg250': {'predicted_species': str(root / 'pg250' / 'predicted_species.parquet'),
                        'primary': {'first': str(root / 'pg250/primary/report-first-pass.parquet'), 'final': str(root / 'pg250/primary/report.parquet'),
                                    'empirical': str(root / 'pg250/primary/report-lib.parquet'), 'design': str(root / 'pg250/design-primary.tsv'),
                                    'exports': str(root / 'pg250/exports-primary.tsv')},
                        'replication': {'first': str(root / 'pg250/replication/report-first-pass.parquet'), 'final': str(root / 'pg250/replication/report.parquet'),
                                        'design': str(root / 'pg250/design-replication.tsv'), 'exports': str(root / 'pg250/exports-replication.tsv')},
                        'within_organism': {sp: {'first': str(root / f'pg250/within-{sp}/report-first-pass.parquet'),
                                                 'final': str(root / f'pg250/within-{sp}/report.parquet'),
                                                 'empirical': str(root / f'pg250/within-{sp}/report-lib.parquet'),
                                                 'design': str(root / f'pg250/design-within-{sp}.tsv'),
                                                 'exports': str(root / f'pg250/exports-within-{sp}.tsv')} for sp in 'HYE'}}}
    import yaml
    path = root / 'config.synthetic.yaml'
    # Relative keys (candidate_counts) must resolve from the repository; write absolute paths for those too.
    config['candidate_counts'] = str(REPO / 'references' / 'candidate-counts.tsv')
    path.write_text(yaml.safe_dump(config, sort_keys=False), encoding='utf-8')
    return path


if __name__ == '__main__':
    print(build(sys.argv[1]))
