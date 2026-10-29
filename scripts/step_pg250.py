"""The 250 pg dataset, one co-search: accepted events, wrong-organism classes and fold changes.

--study primary (all 54 runs, 16 January mixtures): accepted events, retained/added/lost classes of the 18
single-organism runs, A/P counts, two-way bootstrap intervals and contrasts (Figures 1C, 2; Tables S7, D12),
then completeness selection and fold changes (Tables S8, S11; Figure 3). About 25 minutes.
--study replication (9 January mixtures): accepted events and fold changes only (Table S11, D10).
Outputs: <output_dir>/pg250/<study>/export and .../fold-change.
"""
from _setup import setup, note


def main():
    args = setup(__doc__, lambda p: p.add_argument('--study', choices=['primary', 'replication'], required=True))
    from mbraudit import config, pipeline
    block = config.value('pg250')
    roles = ['first', 'final', 'empirical', 'design', 'exports'] if args.study == 'primary' else ['first', 'final', 'design', 'exports']
    job = {role: config.resolve(block[args.study][role]) for role in roles}
    job['predicted_species'] = config.resolve(block['predicted_species'])
    root = config.output_root() / 'pg250' / args.study
    hashes = pipeline.export(job, root / 'export', args.study)
    inputs = {'events': root / 'export' / 'accepted_events.tsv', 'samples': job['design'], 'status': root / 'export' / 'run_status.tsv'}
    pipeline.fold_changes(inputs, root / 'fold-change', args.study)
    note(root / 'step.json', {'step': 'pg250', 'study': args.study, 'mixture_date': pipeline.DATES[args.study], 'input_sha256': hashes,
                              'inputs': {k: str(v) for k, v in job.items()}})
    print('250 pg', args.study, 'complete:', root, flush=True)


if __name__ == '__main__':
    main()
