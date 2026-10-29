# Data

## Raw data

| Dataset | Accession | Article |
| --- | --- | --- |
| Single cells: 64 human cells of four lines (HeLa, HepG2, A549, NCI-H1975), active-matrix digital microfluidic chip, Orbitrap Astral | [PXD048434](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD048434) | Yang et al., JACS Au 2024;4(5):1811-1823, doi:10.1021/jacsau.4c00027 |
| 250 pg digests of HeLa, yeast and E. coli protein extract, alone or mixed, timsTOF Pro 2 diaPASEF | [PXD056832](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD056832) | Wang et al., Nature Communications 2025;16, doi:10.1038/s41467-025-65174-4 |

The single-cell cohort is the first 16 planned cells of each line (`mbraudit/singlecell/design.py`). The 250 pg
analysis uses 54 acquisitions (`design/samples54.tsv`): the 30 three-organism mixtures of 16 January 2025 (six of each
of 4H2Y2E, 5H1Y4E, 5H2Y3E, 5H3Y2E, 5H4Y1E) and, from 12 March 2025, six human-only, six yeast-only and six E. coli-only
runs and six blanks. The 9 January 2025 mixture batch (`design/samples54-replication.tsv`) is used for fold changes only.

## Searches (DIA-NN 2.5.1, two passes)

Every search used DIA-NN 2.5.1 with match-between-runs (`--reanalyse`, empirical-library output), a run precursor
q-value of 0.01 and `--matrices`. Command lines as recorded, with input lists and local paths replaced by braces
(full records in Supplement S1.4, Table S2, and on request):

Single cells, original configuration (refinement on, DIA-NN default):

```
diann.exe --f {64 cell RAW files} --lib {ORIGINAL_PREDICTED_LIBRARY} --threads 27 --qvalue 0.01 --mass-acc 20 --mass-acc-ms1 20 --window 9 --min-pr-charge 2 --max-pr-charge 4 --unimod4 --matrices --reanalyse --out {OUT}/report.tsv --temp {OUT}/tmp --out-lib {OUT}/report-lib.parquet --report-decoys
```

Single cells, refinement off: the same arguments, inputs, input order and library bytes, plus `--no-refine-q`.

Single cells, standard library (refinement off):

```
diann.exe --f {64 cell RAW files} --lib {STANDARD_PREDICTED_LIBRARY} --threads 32 --qvalue 0.01 --mass-acc 20 --mass-acc-ms1 20 --window 9 --min-pr-charge 2 --max-pr-charge 4 --unimod4 --matrices --reanalyse --out {OUT}/report.tsv --temp {OUT}/tmp --out-lib {OUT}/report-lib.parquet --report-decoys --no-refine-q
```

Single-cell pools (Supplement S6): 16 searches of 6 or 16 cells with the original library and the original settings.

250 pg, all 54 runs (automatic mass calibration kept):

```
diann.exe --threads 32 --qvalue 0.01 --f {54 .d folders} --lib {THREE_ORGANISM_LIBRARY} --fasta {HUMAN_YEAST_ECOLI_FASTA} --min-fr-mz 200 --max-fr-mz 1800 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 2 --max-pr-charge 4 --min-pep-len 7 --max-pep-len 30 --cut K*,R* --missed-cleavages 1 --met-excision --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --matrices --reanalyse --no-refine-q --report-decoys --out {OUT}/report.tsv --temp {OUT}/tmp --out-lib {OUT}/report-lib.parquet
```

250 pg, each organism's six single-organism runs on their own (post hoc): the same arguments with the six runs of one
organism and `--use-quant`, reusing the first-pass caches of the 54-run search. 250 pg, 9 January mixtures: the same
arguments with the 9 January mixtures, the same 18 single-organism runs and the six blanks.

### Libraries

| Library | Used by | Generation | Identity |
| --- | --- | --- | --- |
| Original (restricted) predicted library, human plus yeast trap | Single cells, original and refinement off; pools | DIA-NN 2.5.0: trypsin (K/R, not before P), no missed cleavages, no variable modifications, fixed carbamidomethylation, lengths 7 to 30, charges 2 to 4, precursor m/z 200 to 1200, fragment m/z 200 to 1800; 1,745,123 precursors | 423,839,538 bytes, SHA256 7bfbd7c403b8c2e0079e0733c8bc306530086f6a18b486141433bc662f065d53 (a copy whose header path string is shortened, payload byte-identical, was used for the candidate census) |
| Standard predicted library, human plus yeast trap | Single cells, standard library | DIA-NN 2.5.1: one missed cleavage, methionine oxidation and excision, fixed carbamidomethylation, lengths 7 to 30, charges 2 to 4, precursor m/z 300 to 1800, fragment m/z 200 to 1800; 5,878,654 precursors | SHA256 50f32f92449f446b1be2f4e07684c4f42f988a62aa9954845e48023fcada3c28 |
| Three-organism predicted library | 250 pg searches | As the standard library, with human, yeast and E. coli sequences; 6,468,920 precursors | SHA256 8ea00c66e446a505987778382b090593013326f87e9ca88abdfd2d86d12ea38e |
| Empirical libraries | Second pass of each search | Built by DIA-NN from each search's first pass (report-lib.parquet), default inclusion rule | Per search, listed below |

## Processed inputs read by this code (available on request)

The reports and libraries of every search above, the predicted libraries, the FASTA files, the species map and the
verification records are available from the corresponding author on request. The files this repository reads, with the SHA256 prefixes of the
copies used to reproduce the article:

| Configuration key | File | SHA256 (first 16) |
| --- | --- | --- |
| `single_cell.target_fasta` | Human Swiss-Prot FASTA | b205ed9970fee624 |
| `single_cell.entrapment_fasta` | Human plus yeast trap FASTA | e6139c1b385ab288 |
| `single_cell.configurations.on` | report-first-pass.parquet, report.parquet, report-lib.parquet | 5d038f25d457bd84, 4366a343a7fa6816, ce8b9c835693ad98 |
| `single_cell.configurations.off` | the same three files | 8bbdec28022a2759, ae65d6859d8aeb63, 359da591b626d4b3 |
| `single_cell.configurations.standard` | the same three files | 42dc69cc548b8e3a, 2f37560104583813, 47fb05535e0e6d4f |
| `single_cell.pools.*` | 16 pool searches: report-first-pass.parquet, report.parquet | listed in each step.json |
| `libraries.restricted.library` | Original predicted library | see above |
| `pg250.predicted_species` | Species map of the three-organism library (output of the candidate audit) | 03e9fd11871ec8d2 |
| `pg250.primary` | report-first-pass.parquet, report.parquet, report-lib.parquet; report-first-pass.stats.tsv, report.stats.tsv, report.log.txt | 4fe1a2a7f69fee52, d8456303bad9d3a5, be7735bd84a88adf; 08f5c9cc95e29ab4, 501b60a80aa49194, 9f103605ee980223 |
| `pg250.replication` | report-first-pass.parquet, report.parquet | 1c19bb4d4d8dbd12, 69921ded7a3aea1a |
| `pg250.within_organism.H` | report-first-pass.parquet, report.parquet, report-lib.parquet | 22b224a25facc7e3, 9fc4d946ab03e654, dca4b5b934fb8749 |
| `pg250.within_organism.Y` | the same three files | 0a3b0325f4f4927e, 45c212d402507b53, 940bc23516d00273 |
| `pg250.within_organism.E` | the same three files | f594236c60aa5c0f, 92e70ff22a017be3, fbfa6a935fe774ad |

Small inputs are in this repository: the sample sheets and export-status tables of the 250 pg searches (`design/`;
`samples54-replication.tsv` is the 9 January sheet) and the candidate counts of the three predicted libraries
(`references/candidate-counts.tsv`). Every step records the SHA256 of what it read in its `step.json`.
