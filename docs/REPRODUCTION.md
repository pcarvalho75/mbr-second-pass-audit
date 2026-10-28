# Reproduction

## How to reproduce

1. Install Python 3.12 and the packages of `requirements.txt`.
2. Download the processed-data archive (Zenodo, DOI to be added) and copy `config.template.yaml` to
   `config.local.yaml`. Set each path to the archived file (docs/DATA.md lists them with their SHA256 prefixes).
   The configuration file is the only place where local paths appear; it is ignored by git.
3. Run `python scripts/run_all.py`. It runs every step in its own process (about 50 minutes on a desktop with 32 GB of
   memory; the 54-run fold-change step alone takes about 25 minutes, the restricted-library census about 10 minutes;
   `--skip-predicted-census` uses the released candidate counts instead). The 250 pg steps refuse to write into a
   non-empty output folder, so delete the output folder to rerun them.
4. The last step, `scripts/verify.py`, compares every number of the article (`references/published-values.tsv`, 308
   values from the text, Tables 1 to 4, the figure captions and Table S10) with the recomputed value, and writes
   `verification.tsv`. A value matches when it agrees within half a unit of the last printed digit; counts must agree
   exactly. `python scripts/supplementary_data.py --compare supplementary-data-v6.xlsx` checks the published
   Supplementary Data workbook cell by cell.

Each step writes `step.json` with the SHA256 of every input it read. The random draws are fixed (seed 20261018,
2,000 draws), so a rerun gives the same intervals.

The tests (`python -m unittest discover -s tests`, about 2 minutes) need no data: unit tests of the definitions, the
released synthetic fixtures of the fold-change and completeness steps (output hashes must match), hand-counted cases of
the post hoc and within-organism analyses, and a synthetic end-to-end run of every step.

## Verification result

Run on the archived inputs on 28 October 2026 (Windows 11, Python 3.12.10, the versions of `requirements.txt`):

| Check | Result |
| --- | --- |
| Numbers of the text, Tables 1 to 4, figure captions and Table S10 | 308 of 308 match |
| Supplementary Data workbook, 11 sheets | 4,864 of 4,864 compared cells match (the interpretive text columns of D9 and the source-line columns of D9 and D11 are not computed) |
| Intermediate outputs of the original analysis | 350 saved output files compared: 269 byte for byte, 74 identical after line endings, 3 identical after renaming the study labels, 4 with the intended difference listed below |
| Figures | Every plotted value of Figure 1 equals the published figure data; the 2,775 dots of Figure 2 are identical, and its counts, missed-call fractions and medians are among the values checked by `verify.py` |
| Library census | The census of the original predicted library (r = 374,286 / 1,365,729) is identical to the original one |
| Scientific functions | Copied unchanged from the code that produced the published outputs; identity of the syntax trees is recorded for every function (246 definitions) |

### What is not recomputed from the reports

- **Table S1**, the acquisition summary of the single-cell dataset, comes from the instrument file headers, not from
  the search reports. The cohort itself is fixed in `mbraudit/singlecell/design.py`.
- **The candidate ratio of the standard predicted library** (r = 0.2718) and **the species map of the three-organism
  library** are read from the saved output of the library audit (`references/candidate-counts.tsv`; the species map
  in the archive). Regenerating them uses `mbraudit.library.opportunity_audit.audit`, which needs the optional
  pyahocorasick package. The ratio of the original library (r = 0.2741) is recomputed by
  `scripts/step_predicted_census.py` and matches.

### Intended differences from the original code

- Each single-cell configuration now uses its own empirical-library census when its conditional scales are written
  (`conditional_scales.tsv`); the original run reused the refinement-on census for the other two configurations and
  asserted its value. No published number reads that table; the published empirical scale (Table S3, D2) already used
  each configuration's own ratio and is unchanged.
- Input, output and configuration handling, the study labels of the post hoc outputs (`jan16`, `jan09`) and two
  display strings of the within-organism analysis were rewritten; no calculation changed.
