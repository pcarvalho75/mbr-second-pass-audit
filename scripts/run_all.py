"""Run the whole pipeline in order, each step in its own process, then the tables, figures, workbook and check.

    python scripts/run_all.py [--config config.local.yaml] [--skip-predicted-census] [--compare WORKBOOK.xlsx]

Steps (about 50 minutes on a desktop; the 54-run fold-change step alone takes about 25 minutes):
 1. step_single_cell.py --configuration on, off, standard   Tables 1, 2, S3 to S5, D1 to D4; Figure 1
 2. step_compare_configurations.py                          side-by-side tables, empirical-library ratios
 3. step_pools.py                                           pools of one or four cell lines (Table S6, D6)
 4. step_predicted_census.py                                r of the restricted library (D5; optional, slow)
 5. step_pg250.py --study primary, replication              Tables 3, S7, S8, S11, D7, D8, D10; Figure 2
 6. step_within_organism.py --organism H, Y, E              Table 4, Table S9, D9
 7. step_posthoc_fc.py                                      P1 to P3 (Table S11, D11)
 8. step_native_qa.py                                       blank runs (Table S10)
 9. main_tables.py, figure1.py, figure2.py, supplementary_data.py, verify.py
Each step refuses to overwrite a non-empty output folder of the 250 pg steps; delete outputs to rerun.
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--config', type=Path, default=None)
    parser.add_argument('--skip-predicted-census', action='store_true')
    parser.add_argument('--compare', type=Path, default=None, help='published Supplementary Data workbook to check against')
    args = parser.parse_args()
    extra = ['--config', str(args.config.resolve())] if args.config else []
    steps = [['step_single_cell.py', '--configuration', c] for c in ['on', 'off', 'standard']]
    steps += [['step_compare_configurations.py'], ['step_pools.py']]
    if not args.skip_predicted_census:
        steps.append(['step_predicted_census.py'])
    steps += [['step_pg250.py', '--study', 'primary'], ['step_pg250.py', '--study', 'replication']]
    steps += [['step_within_organism.py', '--organism', o] for o in ['H', 'Y', 'E']]
    steps += [['step_posthoc_fc.py'], ['step_native_qa.py'], ['main_tables.py'], ['figure1.py'], ['figure2.py']]
    steps.append(['supplementary_data.py'] + (['--compare', str(args.compare.resolve())] if args.compare else []))
    steps.append(['verify.py'])
    for step in steps:
        print('==>', ' '.join(step), flush=True)
        result = subprocess.run([sys.executable, '-u', str(HERE / step[0]), *step[1:], *extra])
        if result.returncode:
            sys.exit('Step failed: ' + ' '.join(step))
    print('All steps complete.', flush=True)


if __name__ == '__main__':
    main()
