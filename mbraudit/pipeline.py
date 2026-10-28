"""Drivers of the 250 pg steps: the order of calls of the original analysis, with inputs and outputs passed in.

The entry scripts call these with the configured paths; the tests call them with small synthetic fixtures.
"""
from pathlib import Path

import pandas as pd

DATES = {'primary': '20250116', 'replication': '20250109'}


def export(job, out, study, synthetic=False):
    """Accepted events from the reports; for the primary study also classes, A/P counts, bootstrap and contrasts."""
    from mbraudit.pg250 import endpoints as engine
    from mbraudit.pg250.guard import fresh_output, verify_inputs
    inputs, hashes = verify_inputs(job)
    out = fresh_output(out, inputs)
    config = {'inputs': {k: str(v) for k, v in inputs.items()}, 'output': str(out), 'mixture_date': DATES[study]}
    design = pd.read_csv(inputs['design'], sep='\t', keep_default_na=False)
    engine.design_check(design, synthetic, mixture_date=config['mixture_date'])
    events, status, empirical = engine.adapt_reports(config, design, out)
    if study == 'primary':
        members, available = engine.memberships(events, status, design, empirical)
        summary, a, p, pairs, keys = engine.endpoint_summary(members, design, available)
        engine.table(out / 'endpoint_availability.tsv', available)
        engine.table(out / 'fixed_memberships.tsv', members)
        engine.table(out / 'endpoint_counts.tsv', summary)
        engine.table(out / 'bootstrap_event_axis.tsv', pairs)
        engine.table(out / 'bootstrap_key_axis.tsv', pd.DataFrame(keys, columns=engine.KEY))
        if len(summary):
            ledger = engine.bootstrap(summary, a, p, pairs, keys, design, out)
            engine.contrasts(summary, ledger, out)
    return hashes


def fold_changes(inputs, out, study, mode='real', sample_names=()):
    """Completeness selection (primary only) and fold changes from accepted events, samples and run status."""
    from mbraudit.pg250 import fold_change, fold_change_replication
    from mbraudit.pg250.guard import fresh_output, guard_inputs
    engine = fold_change if study == 'primary' else fold_change_replication
    verified = guard_inputs(inputs, mode, producer=engine.__file__, sample_names=sample_names)
    events, samples, status = (pd.read_csv(inputs[k], sep='\t') for k in ['events', 'samples', 'status'])
    if study == 'primary':
        e, s, status, keys = engine.validate(events, samples, status)
    else:
        e, s, status, keys = engine.validate(events, samples, status, mixture_date=DATES[study])
    e = engine.transitions(e)
    tables = {}
    if study == 'primary':
        t, c, m, u = engine.t3(e, s, status, keys, guard=verified)
        tables.update(t3_summary=t, t3_coverage=c, t3_membership=m, t3_unavailable=u)
        f, fs, fu = engine.fc(e, s, status, keys, guard=verified)
    else:
        f, fs, fu = engine.fc(e, s, status, keys, guard=verified, mixture_date=DATES[study])
    keytable = pd.DataFrame(keys, columns=engine.KEY)
    keytable.insert(0, 'key', range(len(keys)))
    tables.update(fc_precursors=f, fc_summary=fs, fc_unavailable=fu, keys=keytable)
    out = fresh_output(out, {k: Path(v).resolve() for k, v in inputs.items()})
    for name, frame in tables.items():
        engine.output_frame(out / (name + '.tsv'), frame)
    return tables


def within_organism(job, out, organism, synthetic=False):
    """One organism's six single-organism runs analysed alone."""
    from mbraudit.pg250 import within_organism as engine
    from mbraudit.pg250.guard import fresh_output, verify_inputs
    inputs, hashes = verify_inputs(job)
    out = fresh_output(out, inputs)
    settings = {'inputs': {k: str(p) for k, p in inputs.items()}, 'input_sha256': hashes, 'output': str(out)}
    design = pd.read_csv(inputs['design'], sep='\t', keep_default_na=False)
    engine.design_check(design, organism, synthetic)
    events, status, empirical = engine._base().adapt_reports(settings, design, out)
    engine.interface_check(events, status, design)
    members, available = engine.memberships(events, status, design, empirical, organism)
    summary, a, p, pairs, keys = engine.endpoint_summary(members, design, available, organism)
    for name, frame in [('endpoint_availability', available), ('fixed_memberships', members),
                        ('endpoint_counts', summary), ('bootstrap_event_axis', pairs),
                        ('bootstrap_key_axis', pd.DataFrame(keys, columns=engine.KEY)),
                        ('planned_endpoint_status', engine.planned_cells(design, available, summary))]:
        engine.table(out / (name + '.tsv'), frame)
    if len(summary):
        ledger = engine.bootstrap(summary, a, p, pairs, keys, design, out)
        engine._base().contrasts(summary, ledger, out)
        engine.interpretation_labels(out)
    return hashes


def posthoc(studies, out, synthetic=False):
    """studies: {alias: {'P3': bool, 'batch_date': 'YYYYMMDD', 'inputs': {role: path}}}, in output order."""
    from mbraudit.pg250 import posthoc_fc as engine
    from mbraudit.pg250.guard import fresh_output, verify_inputs
    scopes, inputs, hashes = {}, {}, {}
    for alias, scope in studies.items():
        names, pins = verify_inputs(scope['inputs'])
        scopes[alias] = {'P3': scope['P3'], 'batch_date': scope['batch_date'],
                         'inputs': {k: {'path': str(p), 'sha256': pins[k]} for k, p in names.items()}}
        inputs.update({alias + '.' + k: p for k, p in names.items()})
        hashes.update({alias + '.' + k: v for k, v in pins.items()})
    out = fresh_output(out, inputs)
    guard = {'seal': engine._SEAL, 'config': {'studies': scopes}, 'inputs': inputs, 'synthetic': synthetic, 'receipt': {'inputs': hashes}}
    tables = engine.analyse(guard)
    for name, frame in tables.items():
        if name == 'source_cells':
            continue  # a cell-by-cell copy of the inputs (about 1 GB for the real data); no result depends on it
        frame.to_csv(out / (name + '.tsv'), sep='\t', index=False, lineterminator='\n', na_rep='NA')
    return hashes
