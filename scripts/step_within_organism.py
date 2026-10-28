"""The 250 pg dataset, post hoc: one organism's six single-organism runs analysed on their own, without the
mixtures (Results 3.4; Table 4; Table S9, D29). Same three-organism database and settings as the 54-run analysis.
Outputs: <output_dir>/pg250/within-<organism>.
"""
from _setup import setup, note


def main():
    args = setup(__doc__, lambda p: p.add_argument('--organism', choices=['H', 'Y', 'E'], required=True))
    from mbraudit import config, pipeline
    from mbraudit.pg250.within_organism import scope
    block = config.value('pg250')
    job = {role: config.resolve(block['within_organism'][args.organism][role]) for role in ['first', 'final', 'empirical', 'design', 'exports']}
    job['predicted_species'] = config.resolve(block['predicted_species'])
    out = config.output_root() / 'pg250' / ('within-' + args.organism)
    hashes = pipeline.within_organism(job, out, args.organism)
    note(out / 'step.json', {'step': 'within-organism', 'organism': args.organism, 'input_sha256': hashes, 'scope': scope(args.organism)})
    print('Within-organism', args.organism, 'complete:', out, flush=True)


if __name__ == '__main__':
    main()
