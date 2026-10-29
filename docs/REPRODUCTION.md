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
4. The last step, `scripts/verify.py`, compares every number of the article (`references/published-values.tsv`, 362
   values: 308 from the text, Tables 1 and 2, Table D12 (Tables 3 and 4 of the earlier manuscript, now the values of
   Figure 2), the figure captions and Table S10, and 54 of the post hoc carryover check
   of Supplement S9) with the recomputed value, and writes `verification.tsv`. A value matches when it agrees within
   half a unit of the last printed digit; counts must agree exactly. `scripts/supplementary_data.py` writes the
   published Supplementary Data workbook itself; with `--compare WORKBOOK.xlsx` it checks a reference workbook cell by
   cell.

Each step writes `step.json` with the SHA256 of every input it read. The random draws are fixed (seed 20261018,
2,000 draws), so a rerun gives the same intervals.

The tests (`python -m unittest discover -s tests`, about 2 minutes) need no data: unit tests of the definitions, the
released synthetic fixtures of the fold-change and completeness steps (output hashes must match), hand-counted cases of
the post hoc and within-organism analyses, and a synthetic end-to-end run of every step.

## Verification result

Run on the archived inputs on 28 October 2026 (Windows 11, Python 3.12.10, the versions of `requirements.txt`):

| Check | Result |
| --- | --- |
| Numbers of the text, Tables 1 and 2, Table D12, figure captions and Table S10 | 308 of 308 match |
| Carryover check of Supplement S9 (post hoc) | 54 of 54 match |
| Supplementary Data workbook, 12 sheets | Written by `scripts/supplementary_data.py`; D12 (the source data of Figure 2) holds Tables 3 and 4 of the earlier manuscript with their printed digits, and D1 to D11 are unchanged by its addition; against the recomputed tables, 10,730 of 10,730 compared cells match (all blocks, including the frozen readings and source lines of D9 and D11); against the earlier published workbook, every cell value is equal; 330 integral counts once shown as 292194.0 are now shown as 292194, the filter rows of D9 and D11 are listed in the order F1, F2, F3, and the README sentence on how the values are produced and the label of the D11 sensitivity block were reworded |
| Intermediate outputs of the original analysis | 350 saved output files compared: 269 byte for byte, 74 identical after line endings, 3 identical after renaming the study labels, 4 with the intended difference listed below |
| Figures | Every plotted value of Figure 1 equals the published figure data; every plotted value of Figure 2 (`figure2-data.tsv`) equals the full-precision rows of Tables 3 and 4 of the earlier manuscript; the 2,775 dots of Figure 3 (the fold-change figure, Figure 2 of the earlier manuscript) are identical, and its counts, missed-call fractions and medians are among the values checked by `verify.py`; Figures S1 and S3 are pixel-identical to the published files and Figure S2 differs only in the filter label of its two panel titles |
| Library census | The census of the original predicted library (r = 374,286 / 1,365,729) is identical to the original one |
| Scientific functions | Copied unchanged from the code that produced the published outputs, except that the filter labels were renamed to match the article (F3->F1, F4->F2, F0->F3; nothing else changed); identity of the syntax trees was recorded for every function (246 definitions) before the renaming, and every output of a full run after the renaming equals that of the run before it, after the same renaming |

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
- The filters are labelled as in the article (F1 primary, F2, F3; docs/METHODS-IN-CODE.md lists the earlier labels),
  the panels of Figure 1 follow the text (completeness B, class-relative fractions C) and the Supplementary Data
  sheets are numbered in the order the Supplement first cites them. Two tables sorted by the label list their rows in
  a new order (`fc_summary.tsv`, `t3_summary.tsv`); no number changed.
- Tables 3 and 4 of the earlier manuscript are now Figure 2 (`scripts/figure2.py`) and Table D12 of Supplementary
  Data 1; `scripts/main_tables.py` writes Tables 1 and 2 only. The fold-change figure is Figure 3 (`scripts/figure3.py`,
  `figure3.*`, `figure3-points.tsv`); its drawing code is unchanged.
