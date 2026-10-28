"""The 250 pg dataset, post hoc: carryover or false matches (Supplement S9). The wrong-organism identifications
added to the six human-only runs (F1, all 54 runs analysed together) against the yeast-only and E. coli-only runs
acquired earlier the same day: run order, organism, intensity relative to the source runs and recurrence, with the
accepted precursors of the blank runs. Reads the two reports of the 54-run search (about 70 MB) and the species
map; uses the native QA of step_native_qa.py for the calibration failures of the blanks when it is present.
Outputs: <output_dir>/pg250/carryover/carryover.tsv (per-run and summary values).
"""
from _setup import setup, note


def main():
    setup(__doc__)
    from mbraudit import config
    from mbraudit.pg250 import carryover
    from mbraudit.pg250.guard import sha
    root = config.output_root() / 'pg250'
    inputs = {'first': config.input_path('pg250.primary.first'), 'final': config.input_path('pg250.primary.final'),
              'predicted_species': config.input_path('pg250.predicted_species'), 'design': config.input_path('pg250.primary.design')}
    native = root / 'native-qa' / 'native-run-qa.tsv'
    out = root / 'carryover'
    table = carryover.run(inputs['first'], inputs['final'], inputs['predicted_species'], inputs['design'], out, native)
    note(out / 'step.json', {'step': 'carryover', 'mask': carryover.MASK, 'input_sha256': {k: sha(p) for k, p in inputs.items()},
                             'native_qa': str(native) if native.exists() else None})
    print('Carryover check:', len(table), 'values;', out / 'carryover.tsv', flush=True)


if __name__ == '__main__':
    main()
