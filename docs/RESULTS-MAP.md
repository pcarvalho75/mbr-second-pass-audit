# Results map

Where each table, figure and number of the article comes from: the command that produces it, the output file (under
the configured `output_dir`) and the value printed in the article. Supplement numbering is that of the final short
Supplement (S1 to S12, Tables S1 to S11, Figures S1 to S3) and Supplementary Data 1 (sheets D1 to D11, numbered in
the order the Supplement first cites them). `scripts/verify.py` checks all 362 values of
`references/published-values.tsv`; `scripts/supplementary_data.py` writes the published workbook itself, and with
`--compare` checks a reference workbook cell by cell.

Steps (run in this order by `scripts/run_all.py`):

| Step | Command | Main outputs |
| --- | --- | --- |
| A | `python scripts/step_single_cell.py --configuration on` (also `off`, `standard`) | `single-cell/<conf>/accounting/*.tsv`, `single-cell/<conf>/intervals/*.tsv`, `single-cell/empirical-census-<conf>.tsv` |
| B | `python scripts/step_compare_configurations.py` | `single-cell/comparison-off/*.tsv`, `single-cell/comparison-standard/*.tsv` |
| C | `python scripts/step_pools.py` | `single-cell/pools/A_*.tsv`, `pool_membership.tsv` |
| D | `python scripts/step_predicted_census.py` | `single-cell/predicted-census.tsv` |
| E | `python scripts/step_pg250.py --study primary` (also `replication`) | `pg250/<study>/export/*.tsv`, `pg250/<study>/fold-change/*.tsv` |
| F | `python scripts/step_within_organism.py --organism H` (also `Y`, `E`) | `pg250/within-<organism>/*.tsv` |
| G | `python scripts/step_posthoc_fc.py` | `pg250/posthoc-fc/summary.tsv`, `membership.tsv` |
| H | `python scripts/step_native_qa.py` | `pg250/native-qa/native-run-qa.tsv` |
| I | `python scripts/step_carryover.py` | `pg250/carryover/carryover.tsv` |
| J | `python scripts/main_tables.py`; `figure1.py`; `figure2.py`; `figures_supplementary.py`; `supplementary_data.py`; `verify.py` | `tables/`, `figures/`, `figures/supplementary/`, `supplementary-data.xlsx`, `supplement-table-sources.xlsx`, `key-numbers.tsv`, `verification.tsv` |

## Main tables

| Table | Needs | Command | Output | Printed values (examples) |
| --- | --- | --- | --- | --- |
| Table 1, accounting of the three single-cell configurations (F1) | A | `scripts/main_tables.py` | `tables/table1.tsv` (source: `accounting.tsv`, `represented_coverage.tsv` of each configuration) | first pass 292,194 / 264,295 / 290,999; final 351,868 / 336,692 / 384,573; retained 250,683; added 101,185; lost 41,511 (31,421 in library, 10,090 outside); represented precursors 14,998 to 11,913; common-row completeness 34.92% to 46.15% (refinement on) |
| Table 2, yeast trap fractions by class | A, B | `scripts/main_tables.py` | `tables/table2.tsv`, `table2-final-result.tsv` (sources: `class_species.tsv`, `class_intervals.tsv`, `class_ratio_intervals.tsv`, `comparison-*/class_species_scales_comparison.tsv`, `references/candidate-counts.tsv`) | retained/on 164/246201, 0.067%, 0.019 to 0.132, 1.000, 0.310; added/on 218/100239, 0.217%, 0.089 to 0.383, 3.265, 1.011; lost/on 142/40915, 0.347%, 0.177 to 0.596, 5.210, 1.613; r = 0.2741 and 0.2718; final result 0.110% and 0.513% (on), 0.140% and 0.649% (off), 0.171% and 0.801% (standard) |
| Table 3, wrong-organism identifications, 54 runs | E | `scripts/main_tables.py` | `tables/table3.tsv` (sources: `pg250/primary/export/endpoint_counts.tsv`, `fraction_intervals.tsv`, two_way) | human retained 79/17,602 = 0.004488 [0.002575, 0.006851]; added 220/6,708 = 0.0328 [0.02537, 0.04179]; lost 22/3,487 = 0.006309 [0.002449, 0.01157]; pooled added 0.1726 [0.1498, 0.1965] |
| Table 4, human-only runs analysed with the other 48 runs and on their own | E, F | `scripts/main_tables.py` | `tables/table4.tsv` (source: `endpoint_counts.tsv` of the 54-run and within-H analyses, `summary.py:accepted_passes`) | first pass 101/21,089 (0.479%) and 108/21,107 (0.512%); final 299/24,310 (1.230%) and 118/25,681 (0.459%); added 220/6,708 and 39/7,700; human identifications 24,011 and 25,563 |

## Figures

| Figure | Needs | Command | Output | Checked values |
| --- | --- | --- | --- | --- |
| Figure 1 (A accounting, B completeness under F1, F2 and F3, C class/retained ratios with two-way intervals) | A, E | `scripts/figure1.py` | `figures/figure1.{svg,pdf,png,tiff}`, every plotted value in `figures/figure1-data.json` | B: F1 own rows 30.44 and 46.15, common rows 34.92 and 46.15; C: 3.265 [1.510, 8.669] and 5.210 [2.013, 19.165] (on), 3.256 [1.575, 8.704] and 6.822 [2.722, 20.931] (off), 2.947 [1.599, 5.782] and 3.432 [1.558, 9.531] (standard), 2.725 [2.244, 3.307] and 0.577 [0.441, 0.737] (250 pg) |
| Figure 2 (fold change against intensity, retained and added ions) | E, G | `scripts/figure2.py` | `figures/figure2.*`, every dot in `figures/figure2-points.tsv`, bin medians in `figure2-binned-medians.tsv` | 709 and 417 yeast, 1,046 and 603 E. coli ions; 14 yeast and 42 E. coli ions in both panels |
| Figures S1 to S3 (S1 added identifications of the 250 pg dataset by support and quantity tercile; S2 entities wholly absent, single cells, F1; S3 completeness filters, 250 pg) | A, E | `scripts/figures_supplementary.py` | `figures/supplementary/FigureS1.{pdf,png,tiff}` to `FigureS3`, every plotted value in `figures/supplementary/figure-s-data.json` | S2: 3,085 precursors, 2,992 peptides, 462 groups wholly absent; S3 human final result 0.0123, 0.005842 and 0.004249 (no filter, 3 of 6, 5 of 6) |

## Numbers in the abstract, Results, Key Points and Conclusions

All are in `key-numbers.tsv` (written by `scripts/verify.py`) with the source file of each.

| Number (paragraph) | Source | Printed value |
| --- | --- | --- |
| Added, lost, first pass, final, retained observations (Abstract, 3.1, Key Points, Conclusions) | Table 1, on | 101,185; 41,511; 292,194; 351,868; 250,683 |
| Distinct peptide ions 14,998 to 11,913; completeness 34.92% to 46.15% (Abstract, 3.1) | `represented_coverage.tsv`, F1 precursor | as printed |
| Lost observations keeping their protein group through another peptide (3.1) | `class_species.tsv`, F1 lost | 34,666 (83.51%) |
| Precursors and human-only protein groups wholly absent (3.1; Table S5) | `entity_species.tsv`, `human_group_losses.tsv` | 3,085 and 432 |
| Yeast 3.3 and 5.2 times more frequent among added and lost (Abstract, 3.2, Key Points) | `class_ratio_intervals.tsv`, two_way | 3.265 [1.510, 8.669]; 5.210 [2.013, 19.165] |
| Final result: lower bound and first-pass-scaled estimate (Abstract, 3.2) | `scales.py:final_single_cell_bound` | 0.110% (382/346,440) and 0.513% |
| Mixed against single-line pools (3.2; Table S6) | `pools/A_composition_uncertainty.tsv` | 1.08 (16 cells), 0.79 (6 cells) |
| 220 of 6,708 added (3.28%), seven times the retained 0.4488% (Abstract, 3.3, Key Points) | Table 3 | as printed; per-run added fractions 2.70% to 4.45% (`endpoint_counts.tsv`, run scope; Table S7) |
| Final human-only wrong-organism fraction (3.3, 3.4) | `summary.py:accepted_passes` | 1.230% (299/24,310) |
| Yeast-only and E. coli-only runs (3.3) | Table 3; first pass from `accepted_passes` | added 29.07% and 24.62%; first pass about 11% and 7% |
| Pooled added/retained and lost/retained (3.3) | `contrast_intervals.tsv`, T1 pooled two_way ratio | 2.725 (2.244 to 3.307); 0.5773 (0.4409 to 0.7367) |
| On their own: 39 of 7,700, 82% fewer (Abstract, 3.4, Key Points) | `within-H/endpoint_counts.tsv`; `summary.py:added_count_comparison` | 39/7,700; 82% |
| At least 5 of 6 runs: 1.230% to 0.4249%, 71.24% kept (3.4; Table S8) | `fold-change/t3_summary.tsv`, `t3_coverage.tsv` | as printed |
| Missed two-fold change (Abstract, 3.5, Key Points) | `posthoc-fc/summary.tsv`, P2 | yeast 49/709 (6.91%) retained, 191/417 (45.80%) added; E. coli 51/1,046 (4.88%), 190/603 (31.51%); wrong sign 0.42% and 9.35% |
| Median yeast log2 fold changes (3.5; Table S11) | `fold-change/fc_summary.tsv` own_context | 1.953, 1.735, 1.570, 1.048; 9 January 2.108, 1.866, 1.620, 1.064 |
| Low-abundance added ions; extra compression (3.5) | `posthoc-fc/summary.tsv`, P3 low tercile | 304 of 417 yeast, 434 of 603 E. coli; 0.3374 and 0.3007 |
| Added values with the smaller error in three of four comparisons; MAE lowered by 0.0620 and 0.0563 (3.5) | `fc_summary.tsv` R_T_intersection (both batches); `summary.tsv` P1 | 3 of 4; -0.0620, -0.0563 |
| Carryover or false matches, human-only runs of the 54-run search (post hoc; main text and Supplement S9) | `pg250/carryover/carryover.tsv` (step I) | yeast additions 11, 23, 19, 17, 16, 13 and E. coli additions 17, 28, 24, 14, 14, 24 in acquisition order (12743 to 12749), of 1,032 to 1,187 added; added wrong-organism fractions 2.71% to 4.45%; 121 E. coli against 99 yeast; median ratio to the source runs 0.1410 and 0.1824 (log10 interquartile ranges 1.32 and 1.06; Pearson r 0.06 and 0.16; 88 of 99 and 121 of 121 accepted in the source runs); 167 distinct precursors, 135, 20, 7, 3 and 2 added in 1, 2, 3, 4 and 6 runs; blank 12740: no first-pass identification, 4 calibration failures, 148 final (117 human, 14 yeast, 12 E. coli, 5 shared) |
| Support contrasts (Discussion; Supplement S3) | `intervals/support_contrasts.tsv` (single cells); `export/contrasts.tsv` T2 (250 pg) | 5.798 (105/17,423 against 22/21,165); 2.95 (214/593 against 304/2,485) |
| 54 acquisitions: 30 mixtures, 18 single-organism runs, 6 blanks (Methods 2.4) | `design/samples54.tsv` | as printed |

## Short Supplement tables

| Table | Content | Source (Supplementary Data sheet) |
| --- | --- | --- |
| Table S1 | Cell acquisitions and cohort membership | Acquisition manifest (RAW headers); not recomputed here. The cohort itself is `mbraudit/singlecell/design.py`. |
| Table S2 | The searches and command lines | docs/DATA.md (recorded command lines) |
| Table S3 | Observation classes of the three configurations with scales | `supplement-table-sources.xlsx` sheet "Table S3 (v5.1 S72)"; D3, D4 |
| Table S4 | Accounting and completeness, F1/F2/F3, on and off | D2, D5 |
| Table S5 | Entities wholly absent | `supplement-table-sources.xlsx` sheet "Table S5 (v5.1 S71)" |
| Table S6 | Mixed/homogeneous pools | `supplement-table-sources.xlsx` sheet "Table S6 (v5.1 S2)"; D6 |
| Table S7 | Human-only runs, per run | D7 |
| Table S8 | Completeness filters, 250 pg | D8 |
| Table S9 | Within-organism co-searches | D9 |
| Table S10 | Blank runs, native QA | `supplement-table-sources.xlsx` sheet "Table S10 (v5.1 S56)" (`native-run-qa.tsv`) |
| Table S11 | Fold changes by context and batch | D10, D11 |

## Supplementary Data 1 (`scripts/supplementary_data.py`)

| Sheet | Content | Source outputs |
| --- | --- | --- |
| D1 | Candidate ratios: predicted and empirical libraries (A) and empirical census on/off (B) | `predicted-census.tsv`, `empirical-census-on.tsv`, `comparison-off/empirical_library_opportunities.tsv` |
| D2 | Accounting of the three configurations, F1/F2/F3 | `single-cell/<conf>/accounting/accounting.tsv` |
| D3 | Class composition, ratios, scales, group representation | `single-cell/comparison-off`, `comparison-standard/class_species_scales_comparison.tsv` |
| D4 | Class-fraction (A) and class-ratio (B) intervals | `intervals/class_intervals.tsv`, `class_ratio_intervals.tsv` |
| D5 | Represented entities and completeness | `accounting/represented_coverage.tsv` |
| D6 | Pool definitions (A) and per-pool results (B) | `pools/pool_membership.tsv`, `pools/A_per_pool.tsv` |
| D7 | Per-run accounting, 18 single-organism runs | `pg250/primary/export/endpoint_counts.tsv`, `fraction_intervals.tsv` |
| D8 | Completeness filters | `pg250/primary/fold-change/t3_summary.tsv`, `t3_coverage.tsv` |
| D9 | Within-organism contrasts (A, with the frozen reading rule and the source line of each interval) and class fractions (B) | `pg250/within-*/contrasts.tsv`, `contrast_intervals.tsv`, `endpoint_counts.tsv` |
| D10 | Fold changes by batch (A, B), 16 January contexts (C) and intersections (D) | `pg250/<study>/fold-change/fc_summary.tsv` |
| D11 | P1 to P3 with sensitivities (source line in `summary.tsv`) | `pg250/posthoc-fc/summary.tsv` |
