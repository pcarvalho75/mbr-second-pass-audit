"""Copy DIA-NN's native per-run statistics (report*.stats.tsv) onto the designed runs.

Restates the old analysis/native_run_qa.py main() as a function (canonical_sample is unchanged). Used for the
blank-run QA counts (Supplement S9, Table S10). The realised calibration tolerances and calibration failures
of that table are read from DIA-NN's log by log_calibration().
"""
from pathlib import Path, PurePosixPath

import pandas as pd

COLUMN_MAP = {
    'Precursors.Identified': 'Precursors.Identified',
    'Median.Mass.Acc.MS1': 'Median.Mass.Acc.MS1',
    'Median.Mass.Acc.MS1.Corrected': 'Median.Mass.Acc.MS1.Corrected',
    'Median.Mass.Acc.MS2': 'Median.Mass.Acc.MS2',
    'Median.Mass.Acc.MS2.Corrected': 'Median.Mass.Acc.MS2.Corrected',
}


def canonical_sample(value):
    return PurePosixPath(str(value).replace('\\', '/')).stem


def native_qa(design_path, sources, column_map=COLUMN_MAP):
    """sources: list of {'pass': 'first'|'final', 'path': stats.tsv}. Returns one row per designed run and pass."""
    design = pd.read_csv(design_path, sep='\t', keep_default_na=False)
    if design['sample'].duplicated().any():
        raise ValueError('Designed sample names must be unique')
    tables = []
    for spec in sources:
        frame = pd.read_csv(Path(spec['path']), sep='\t', keep_default_na=False)
        sample_column = spec.get('sample_column', 'File.Name')
        frame['sample'] = frame[sample_column].map(canonical_sample)
        if frame['sample'].duplicated().any() or not set(frame['sample']) <= set(design['sample']):
            raise ValueError('Native stats contain duplicate or undesigned runs')
        selected = frame[['sample']].copy()
        for output, source in column_map.items():
            if source not in frame:
                raise ValueError('Native statistic is unavailable: ' + source)
            selected[output] = frame[source]
        merged = design.merge(selected, on='sample', how='left', validate='one_to_one')
        merged.insert(0, 'pass', spec['pass'])
        merged.insert(1, 'manifest_ordinal', range(1, len(merged) + 1))
        tables.append(merged)
    return pd.concat(tables, ignore_index=True)


def log_calibration(log_path):
    """Per pass and run, the realised calibration tolerances and the number of failed calibration attempts,
    read from DIA-NN's log (lines 'File #k/n', 'Loading run', 'Calibrating with mass accuracies A (MS1), B (MS2)',
    'Cannot perform ... mass calibration'). The second pass starts at the line 'Second pass:'.
    """
    import re
    rows, current, stage = [], None, 'first'
    for line in Path(log_path).read_text(encoding='utf-8', errors='replace').splitlines():
        if line.startswith('Second pass:'):
            stage = 'final'
            continue
        m = re.match(r'^\[[0-9:]+\] File #(\d+)/(\d+)', line)
        if m:
            current = {'pass': stage, 'manifest_ordinal': int(m.group(1)), 'sample': None,
                       'search_ms1_ppm': None, 'search_ms2_ppm': None, 'calibration_failures': 0}
            rows.append(current)
            continue
        if current is None:
            continue
        m = re.match(r'^\[[0-9:]+\] Loading run (.+)$', line)
        if m:
            current['sample'] = canonical_sample(m.group(1).strip())
            continue
        m = re.search(r'Calibrating with mass accuracies ([0-9.]+) \(MS1\), ([0-9.]+) \(MS2\)', line)
        if m:
            current['search_ms1_ppm'], current['search_ms2_ppm'] = float(m.group(1)), float(m.group(2))
            continue
        if re.search(r'Cannot perform (MS1 )?mass calibration', line):
            current['calibration_failures'] += 1
    return pd.DataFrame(rows)
