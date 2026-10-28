"""The 250 pg dataset: DIA-NN's native per-run statistics of the 54-run search, both passes
(blank-run QA: identified precursors, realised calibration tolerances and calibration failures; Supplement S9,
Table S10). Statistics come from report*.stats.tsv, tolerances and failures from report.log.txt.
Outputs: <output_dir>/pg250/native-qa/native-run-qa.tsv.
"""
from _setup import setup, note


def main():
    setup(__doc__)
    from mbraudit import config
    from mbraudit.pg250.guard import sha
    from mbraudit.pg250.native_qa import native_qa, log_calibration

    block = config.value('pg250.primary')
    design = config.resolve(block['design'])
    sources = [{'pass': 'first', 'path': config.resolve(block['stats_first'])},
               {'pass': 'final', 'path': config.resolve(block['stats_final'])}]
    out = config.output_dir('pg250', 'native-qa')
    log = config.resolve(block['log'])
    qa = native_qa(design, sources).merge(log_calibration(log), on=['pass', 'manifest_ordinal', 'sample'], how='left', validate='one_to_one')
    qa.to_csv(out / 'native-run-qa.tsv', sep='\t', index=False, na_rep='NA', lineterminator='\n')
    note(out / 'step.json', {'step': 'native-qa', 'input_sha256': {**{s['pass']: sha(s['path']) for s in sources}, 'log': sha(log)},
                             'design_sha256': sha(design)})
    print('Native QA complete:', out, flush=True)


if __name__ == '__main__':
    main()
