"""Single-cell dataset: species census of the restricted predicted library, the search space of the original and
refinement-off searches. Gives r = yeast-only / human-only precursor candidates under the one-way rule (Methods
2.3; Table 2 caption, r = 0.2741) and the peptide-unit census (the precursor rows are Supplementary Data D1).
Decodes the 424 MB binary library and labels about 1.7 million precursors; expect 10 to 20 minutes.
Outputs: <output_dir>/single-cell/predicted-census.tsv.
"""
from _setup import setup


def main():
    setup(__doc__)
    from mbraudit import config
    from mbraudit.singlecell import candidates
    from mbraudit.singlecell.census import census_rows
    from mbraudit.species import SpeciesRegistry
    out = config.output_dir('single-cell', 'predicted-census')
    (out / 'cache').mkdir(exist_ok=True)
    candidates.PRED = config.input_path('libraries.restricted.library')
    candidates.ROOT = out
    registry = SpeciesRegistry(str(config.input_path('single_cell.target_fasta')), str(config.input_path('single_cell.entrapment_fasta')))
    frame = candidates.label(candidates.predicted(), registry, 'predicted')
    census_rows(frame, 'predicted').to_csv(config.output_root() / 'single-cell' / 'predicted-census.tsv', sep='\t', index=False, lineterminator='\n')
    print('Predicted-library census complete', flush=True)


if __name__ == '__main__':
    main()
