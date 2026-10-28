# Methods in code

Each definition of the Methods (section 2) and of Supplement S2 (Definitions), with the function that implements it.
Paths are relative to the repository root; line numbers refer to this version. These functions are copied unchanged
from the code that produced the published numbers (build record: AST identity with the source, see
docs/REPRODUCTION.md), except that the filter labels were renamed to match the article (F3->F1, F4->F2, F0->F3;
nothing else changed; see Filter labels below).

## Reports and acceptance

| Definition | Where |
| --- | --- |
| A report row is a target (Decoy = 0); run names lose their folder, extension and the dataset prefix. | `mbraudit/species.py:accept_rows` (105), `run_name` (31); 250 pg: `mbraudit/pg250/endpoints.py:adapt_reports` (86), `canonical_sample` (50) |
| Observation: one precursor (modified sequence plus charge) in one cell or run. | `mbraudit/singlecell/reports.py:Reports.load` (61, key = (Modified.Sequence, Precursor.Charge)); `mbraudit/pg250/endpoints.py:memberships` (189) |
| F1 (primary): run, global and library precursor q-values all at most 0.01. F2: F1 plus the run, global and library protein-group q-values at most 0.01. F3: run precursor q-value at most 0.01. Missing fields fail. The same rule is applied to each pass separately. | `mbraudit/filters.py:filter_mask` (12; missing values are filled with 2 before the call, `accounting.py:filtered_pair` 212); 250 pg: `endpoints.py:MASKS` and `adapt_reports` |

## Classes (Methods 2.2)

| Definition | Where |
| --- | --- |
| Retained: accepted in both passes; added (second-pass-only): final result only; lost (first-pass-only): first pass only. First pass = retained + lost; final = retained + added (asserted). | `mbraudit/singlecell/accounting.py:all_cell_accounting` (283); 250 pg: `endpoints.py:memberships` (classes D, F, R, T, L) |
| Lost split: in the empirical library (target entries of report-lib.parquet) or outside it. | `accounting.py:all_cell_accounting`; `endpoints.py:memberships` (`library_state`) |
| Represented precursors (rows) and completeness: own-row = positive finite cell-entity entries / (each pass's rows x cells); common-row = the same over the rows present in both passes. Peptide and protein-group levels use a fixed mapping from the union of accepted precursors. | `accounting.py:coverage` (56), `maps_from_union` (35); `reports.py:annotations`, `group_units` |
| Lost observations whose protein group stays in the same cell through another peptide. | `accounting.py:group_representation` (224) |
| Entities wholly absent from the final result, and human-only groups among them (Table S5). | `accounting.py:entity_composition` (243), `fixed_species` (171) |

## Organism calls and error estimates (Methods 2.3, 2.4)

| Definition | Where |
| --- | --- |
| Single cells, one-way rule: human-only accessions = target; yeast-only = entrapment; mixed = shared; unknown accessions = unknown. A yeast-only peptide found in a human protein with I read as L becomes I/L-shared; human peptides are not checked against yeast. | `mbraudit/species.py:classify` (58), `SpeciesRegistry.labels` (81) |
| Observed yeast fraction E/(T+E); shared, I/L-shared and unknown excluded from the denominator. | `accounting.py:species_record` (149) |
| Ratio of class fractions (r-free). | `accounting.py:fraction_ratio` (157); intervals `bootstrap.py:ratio` (19) |
| Combined estimate E(1+1/r)/(T+E), an upper estimate only for identifications accepted in the first-pass search space, with r = yeast / human precursor candidates of the predicted library searched (one-way rule): 374,286 / 1,365,729 (original library), 1,254,857 / 4,617,329 (standard). Whole final result: retained + added counts summed first. | `mbraudit/singlecell/scales.py:scale_point` (9), `final_single_cell_bound` (24); candidate counts: `mbraudit/singlecell/candidates.py:predicted`, `label` with `census.py:census_rows` (restricted library), or `mbraudit/library/opportunity_audit.py:audit` (156; all libraries, needs pyahocorasick) |
| Historical r = 0.276123 and the empirical-library r of each configuration (Supplement S2, Table S3). | `species.py:ENTRAP_RATIO`; `configurations.py:empirical_opportunities` (17), `comparison_tables` (41) |
| Two-way intervals: 2,000 draws (seed 20261018) resampling cells and precursor blocks with multiplicities; one draw shared by all classes; percentile 2.5 and 97.5 with linear interpolation; any undefined draw makes the interval unavailable. | `mbraudit/singlecell/bootstrap.py:multiplicities` (14), `interval` (35); `intervals.py:joint_draws` (32), `class_intervals` (53) |
| 250 pg: A counts identifications exclusive to an organism absent from the run, P those exclusive to the intended organism, fraction A/n with n = A + P; shared, I/L-shared, unknown and quarantined excluded. Organism calls come from the fixed species map of the three-organism predicted library (symmetric exact and I/L checks). | `mbraudit/pg250/endpoints.py:endpoint_summary` (270); species map: `library/opportunity_audit.py:audit`, `primary_label` (100) |
| 250 pg intervals: two-way run (stratified by organism, six runs each) and precursor resampling, 2,000 shared draws; contrasts (added/retained, lost/retained, high/low support) from the same draws. | `endpoints.py:bootstrap` (327), `contrasts` (363) |
| Within-organism co-searches: the same rules with six runs and one stratum. | `mbraudit/pg250/within_organism.py:memberships` (103), `endpoint_summary` (113), `bootstrap` (123) |
| Whole passes of the human-only runs (Table 4): first pass = R + L counts, final = R + T counts. | `mbraudit/pg250/summary.py:accepted_passes` (7), `added_count_comparison` (33) |

## Support and quantity terciles (Supplement S2, S3)

| Definition | Where |
| --- | --- |
| Support of an added observation: the other cells (runs) with a first-pass acceptance of the same precursor, k/(n-1). Single cells: zero, (0,0.2], (0.2,0.5], (0.5,1]; original bins (0,0.2], (0.2,0.4], (0.4,0.6], (0.6,1]. | `accounting.py:filtered_pair` (212), `support_bin` (331); `support.py:original_bin` (16) |
| Single-cell quantity terciles: cut points from each cell's own first-pass positive quantities; an added observation's final quantity is placed against them. | `accounting.py:quantity_terciles` (338) |
| Single-cell support contrast (low against high support in the high tercile, 5.798) with Wilson intervals. | `support.py:support_intervals` (24), `wilson` (7) |
| 250 pg support (other 53 runs; bins 0, (0,0.2], ..., (0.6,1]) and terciles of added observations within each run. | `endpoints.py:support_bin` (57), `quantity_terciles` (67) |

## Completeness filters and fold changes (Methods 2.4, Results 3.4, 3.5)

| Definition | Where |
| --- | --- |
| Completeness filter: keep a precursor present (positive quantity) in at least 3 or 5 of an organism's six runs; fractions and kept shares (Table S8). | `mbraudit/pg250/fold_change.py:t3` (85) |
| Fold change of a peptide ion: median log2 quantity in the six 5H4Y1E runs minus that in the six 5H1Y4E runs, with at least two positive quantities in each; expected +2 for yeast, -2 for E. coli; contexts first pass (D), final (F), retained (R), added (T). | `fold_change.py:fc` (193); 9 January batch: `fold_change_replication.py:fc` |
| Missed two-fold call: absolute log2 fold change below 1; wrong sign: opposite sign to the design; P1 paired retained/final keys; P3 abundance terciles fixed once on the retained/added union from the median log2 final Precursor.Normalised over the twelve runs. | `mbraudit/pg250/posthoc_fc.py:analyse` (37) |
| Figure 2 coordinates: y = class fold change, x = the P3 abundance divided by log2(10). | `mbraudit/figures/figure2.py:points` |

## Carryover or false matches (Supplement S9, post hoc)

| Definition | Where |
| --- | --- |
| Accepted observations of the two reports of the 54-run search under F1 (targets, run, global and library precursor q-values at most 0.01), with the organism call of the fixed species map, classed retained or added per run. | `mbraudit/pg250/carryover.py:accepted` (31), with `filters.py:filter_mask` and `endpoints.py:memberships` |
| Run order and organism: added yeast and E. coli identifications per human-only run in acquisition order (design column acquisition_id). Intensity: each added wrong-organism ion's Precursor.Quantity over its median Precursor.Quantity in the same-day single-organism runs of its organism (final result), summarised by the median (10 to the median log10 ratio) and the interquartile range of the log10 ratios, with the Pearson correlation of the log10 intensities. Recurrence: the number of human-only runs in which each wrong-organism precursor is added. Blanks: accepted target precursors by organism in each pass, with the calibration failures of the native QA. | `carryover.py:analyse` (51) |

## Pools of one or four cell lines (Supplement S6)

| Definition | Where |
| --- | --- |
| Fixed pool memberships. | `mbraudit/singlecell/design.py` |
| Ratio of yeast per added observation, mixed over homogeneous pools; count-model interval (Clopper-Pearson on the yeast split) and cell, precursor-block and two-way bootstrap intervals. | `mbraudit/singlecell/composition.py:composition` (43), `rate_ratio_ci` (89); `bootstrap.py:bootstrap` (77, kind 'composition') |

## Filter labels

The filters are numbered in the order the article introduces them. The code that produced the published numbers,
and outputs written by earlier versions of this repository, used other labels; the masks are the same.

| Label in the article and this version | Label in earlier outputs | Precursor accepted when |
| --- | --- | --- |
| F1 (primary) | F3 | run, global and library precursor q-values all at most 0.01 |
| F2 | F4 | F1 and the run, global and library protein-group q-values at most 0.01 |
| F3 | F0 | run precursor q-value at most 0.01 |

The renaming also applies to labels built from them: workflow names (all_cells_F3 was all_cells_F0), column names
(second_pass_only_with_F3_firstpass was second_pass_only_with_F0_firstpass) and the file names of the workflow
draws. Two tables are sorted by the label and therefore list their rows in a new order with the same content:
fold-change/fc_summary.tsv and fold-change/t3_summary.tsv.
