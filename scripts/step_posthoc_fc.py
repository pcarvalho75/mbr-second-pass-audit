"""The 250 pg dataset, post hoc fold-change analyses (Results 3.5; Table S11; Supplementary Data D34):
P1 paired retained/final keys, P2 call consequences per class (missed two-fold, wrong sign, within 0.5),
P3 abundance terciles fixed on the retained/added union (16 January only).
Needs step_pg250.py for both studies. Outputs: <output_dir>/pg250/posthoc-fc.
"""
from _setup import setup, note

# Study aliases become the 'batch' column of the outputs; order is that of the original analysis.
STUDIES = {'jan16': ('primary', True, '20250116'), 'jan09': ('replication', False, '20250109')}


def main():
    setup(__doc__)
    from mbraudit import config, pipeline
    root = config.output_root() / 'pg250'
    studies = {}
    for alias, (study, p3, date) in STUDIES.items():
        inputs = {'fc': root / study / 'fold-change' / 'fc_precursors.tsv', 'keys': root / study / 'fold-change' / 'keys.tsv'}
        if p3:
            inputs['events'] = root / study / 'export' / 'accepted_events.tsv'
            inputs['samples'] = config.resolve(config.value('pg250.' + study + '.design'))
        studies[alias] = {'P3': p3, 'batch_date': date, 'inputs': inputs}
    hashes = pipeline.posthoc(studies, root / 'posthoc-fc')
    note(root / 'posthoc-fc' / 'step.json', {'step': 'posthoc-fc', 'input_sha256': hashes, 'studies': {k: v[0] for k, v in STUDIES.items()}})
    print('Post hoc fold changes complete:', root / 'posthoc-fc', flush=True)


if __name__ == '__main__':
    main()
