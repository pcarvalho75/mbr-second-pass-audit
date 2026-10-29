# Completeness at a cost: auditing the second pass of match-between-runs in single-cell DIA proteomics

Analysis code for the article of the same title (Briefings in Bioinformatics, Case study; Martins, Fischer, dos Santos,
Sales, Luz, Muotri and Carvalho). It recomputes every table, figure and number of the article from the DIA-NN reports and the other processed inputs.

## The question

Single cells hold so little material that a peptide seen in one cell is often missed in the next. Match-between-runs
fills these gaps. In DIA-NN it works in two passes: a first pass searches each run on its own, a library is built from
the peptides found in any run, and a second pass searches every run again against that library. The second pass is the
final result. It adds identifications, but it can also add peptides to runs that never contained them, and it can drop
some first-pass identifications.

The article asks, identification by identification: what does the second pass add, what does it lose, and how many of
those changes are wrong? It answers with two public datasets in which part of the truth is known, so that errors can
be counted rather than assumed.

## The two datasets

| Name in the article | Accession | Source | What makes errors visible |
| --- | --- | --- | --- |
| The single-cell dataset | [PXD048434](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD048434) | Yang et al., AM-DMF-SCP, JACS Au 2024, [doi:10.1021/jacsau.4c00027](https://doi.org/10.1021/jacsau.4c00027) | 64 human cells (HeLa, HepG2, A549, NCI-H1975; Orbitrap Astral) searched against human plus yeast trap sequences: every accepted yeast peptide is false. |
| The 250 pg dataset | [PXD056832](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD056832) | Wang et al., Nature Communications 2025, [doi:10.1038/s41467-025-65174-4](https://doi.org/10.1038/s41467-025-65174-4) | 250 pg of HeLa, yeast and E. coli protein extract, alone or mixed in known ratios (timsTOF Pro 2, diaPASEF), 54 runs analysed together: any yeast or E. coli identification in a human-only run is a wrong-organism match, and the known ratios test quantification. |

All searches used DIA-NN 2.5.1 in two-pass mode. The settings, command lines and libraries are in
[docs/DATA.md](docs/DATA.md).

## What the code does

1. **Read both reports of each search.** The first-pass report and the final report, with targets only.
2. **Apply the same filter to both passes.** F1, the primary filter, keeps a precursor when its run, global and library
   q-values are all at most 0.01. F2 adds the three protein-group q-values; F3 uses the run q-value alone.
3. **Sort every observation into a class.** An observation is one precursor (modified sequence plus charge) in one cell
   or run. Retained: accepted in both passes. Added: final result only. Lost: first pass only, split by whether it is in
   the empirical library.
4. **Count represented precursors and completeness** (filled entries over cells times rows).
5. **Call the organism of every observation** from its protein annotation (single cells: human, yeast, shared or
   unknown, with a yeast peptide found in a human protein counted as shared; 250 pg: human, yeast or E. coli).
6. **Estimate errors by class.** Single cells: the yeast fraction of each class, its ratio to the retained class, two-way
   (cell and precursor) bootstrap intervals from 2,000 draws, and a combined estimate that scales the yeast
   fraction by the yeast-to-human ratio of the search space (an upper estimate only for identifications accepted in
   the first-pass search space; the second pass searches the empirical library, so for added observations and the final
   result only the lower bound holds). 250 pg: the wrong-organism fraction of each class in the
   single-organism runs, with two-way intervals.
7. **Compare analyses of different sets of runs.** The six human-only runs analysed with the other 48 runs and on
   their own; cell-line pools of one line against pools of four lines; and, post hoc, whether the wrong-organism
   identifications added to the human-only runs look like carryover from the runs acquired before them.
8. **Test fold changes.** Yeast and E. coli differ four-fold between two of the mixtures; the code computes each
   peptide ion's fold change by class and how often it misses a two-fold change.
9. **Write the tables, the main and supplementary figures and Supplementary Data**, and check every number against
   the article.

## Quick start

Requirements: Python 3.12 with the packages in [requirements.txt](requirements.txt). The optional package
pyahocorasick is needed only to regenerate the candidate census of the predicted libraries with the published audit
step; the released census is in references/candidate-counts.tsv.

```
python -m pip install -r requirements.txt
python -m unittest discover -s tests          # unit tests and a small synthetic end-to-end run (about 2 minutes)
```

To reproduce the article from the processed inputs (search outputs and libraries, available from the corresponding author on request), copy
`config.template.yaml` to `config.local.yaml`, set the paths to the inputs ([docs/DATA.md](docs/DATA.md)), then run

```
python scripts/run_all.py        # about 50 minutes; outputs go to the configured output_dir
```

`run_all.py` runs each step in turn and ends with `verify.py`, which compares 362 numbers with the article and reports
any difference. Each step can also be run alone; [docs/RESULTS-MAP.md](docs/RESULTS-MAP.md) lists the step, output
file and expected value of every table, figure and number.

## Layout

| Path | Content |
| --- | --- |
| `mbraudit/` | The package. `species.py`, `filters.py`: acceptance, organism calls and the F1/F2/F3 filters. `singlecell/`: reports, classes and completeness (`accounting.py`), bootstrap intervals (`bootstrap.py`, `intervals.py`), support (`support.py`), configurations side by side (`configurations.py`), pools (`composition.py`), upper bounds (`scales.py`), library census (`candidates.py`, `census.py`). `pg250/`: the 250 pg dataset (`endpoints.py`, `fold_change*.py`, `within_organism.py`, `posthoc_fc.py`, `summary.py`, `native_qa.py`, `carryover.py`). `library/`: readers of DIA-NN predicted libraries. `figures/`, `tables.py`, `workbook.py`, `key_numbers.py`: outputs. |
| `scripts/` | One entry script per analysis step, table set or figure, and `run_all.py`. |
| `design/` | Sample sheets of the 250 pg searches and their export-status tables. |
| `references/` | The candidate counts of the predicted libraries and the published values checked by `verify.py`. |
| `tests/` | Unit tests, released synthetic fixtures and a synthetic end-to-end run. |
| `docs/` | Results map, methods in code, data and reproduction notes. |

The scientific functions are copied unchanged from the code that produced the published numbers, except that the
filter labels were renamed to match the article (F3->F1, F4->F2, F0->F3; nothing else changed); only input, output
and configuration handling was rewritten. [docs/METHODS-IN-CODE.md](docs/METHODS-IN-CODE.md) maps each definition of
the Methods to its function, and [docs/REPRODUCTION.md](docs/REPRODUCTION.md) records the reproduction result.

## Citation and licence

Cite the article (details to be added on publication; see [CITATION.cff](CITATION.cff)) and, when reusing the raw data,
the two dataset papers above. The code is released under the MIT licence ([LICENSE](LICENSE)).
