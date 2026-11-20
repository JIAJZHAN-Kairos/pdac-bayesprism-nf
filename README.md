# BayesPrism bulk RNA deconvolution

Nextflow workflow that runs [BayesPrism](https://github.com/Danko-Lab/BayesPrism) 2.2.3
(commit `19052052`) once per bulk RNA source against a single-cell reference.
The container is the Bioconductor 3.21 image (R 4.5); BayesPrism and its
dependencies are installed into a task-local library at the start of each task
(see `install.log`).

| Parameter | Input |
| --- | --- |
| `mixtures` | Glob of RDS files, each a samples x genes raw count matrix (`mixture_<source>.rds`). |
| `reference` | RDS list with `sc.dat` (cells x genes raw counts), `cell.type.labels`, `cell.state.labels` and `key` (`"Malignant"`). |
| `outdir` | Output directory or S3 prefix. |

Each source is run with `input.type = "count.matrix"`, `key = "Malignant"`,
`outlier.cut = 0.01`, `outlier.fraction = 0.1`. Outputs per source: initial
cell-state and cell-type fractions, final (updated) cell-type fractions,
posterior CVs, the full BayesPrism object and a run receipt.
Sources are deconvolved separately so the malignant reference is never updated
across quantification pipelines.
