"""Species census of a library, in the layout of the original candidate-ratio table.

This restates, for the empirical library only, the counting loop of the old analysis/candidate_ratios.py main()
(the loop was inline there, so there is no function to carry). Labels come from the carried label() function.
Output columns: library, unit, charge, classification, total_candidates, human_only, yeast_only,
shared_accessions, shared_IL, unknown, mixed_annotation, r (yeast_only / human_only).
"""
import numpy as np
import pandas as pd

from mbraudit.singlecell.candidates import label


def empirical_census(registry, library_path, role='empirical'):
    raw = pd.read_parquet(library_path, columns=['Modified.Sequence', 'Precursor.Charge', 'Stripped.Sequence', 'Protein.Ids', 'Decoy'])
    raw = raw[raw.Decoy == 0].drop(columns='Decoy').drop_duplicates()
    assert not raw.duplicated(['Modified.Sequence', 'Precursor.Charge']).any()
    return census_rows(label(raw, registry, role), role)


def census_rows(df, role):
    """Counts by unit (precursor, stripped peptide), charge and classification of a labelled candidate table."""
    output = []
    for unit, col in [('precursor', 'Modified.Sequence'), ('stripped_peptide', 'Stripped.Sequence')]:
        for charge in ['all'] + sorted(df['Precursor.Charge'].unique().tolist()):
            block = df if charge == 'all' else df[df['Precursor.Charge'] == charge]
            if unit != 'precursor':
                block = block.groupby(col, sort=False).agg(species=('species', lambda x: next(iter(set(x))) if len(set(x)) == 1 else 'mixed_annotation'), species_before_IL=('species_before_IL', lambda x: next(iter(set(x))) if len(set(x)) == 1 else 'mixed_annotation'))
            for stage in ['species_before_IL', 'species']:
                counts = block[stage].value_counts()
                target, entrap = (int(counts.get('target', 0)), int(counts.get('entrapment', 0)))
                row = {'library': role, 'unit': unit, 'charge': charge, 'classification': 'accession_only' if stage == 'species_before_IL' else 'I_as_L_exclusive', 'total_candidates': len(block), 'human_only': target, 'yeast_only': entrap, 'shared_accessions': int(counts.get('shared', 0)), 'shared_IL': int(counts.get('shared_IL', 0)), 'unknown': int(counts.get('unknown', 0)), 'mixed_annotation': int(counts.get('mixed_annotation', 0)), 'r': entrap / target if target else np.nan}
                assert sum((row[k] for k in ['human_only', 'yeast_only', 'shared_accessions', 'shared_IL', 'unknown', 'mixed_annotation'])) == len(block)
                output.append(row)
    return pd.DataFrame(output)
