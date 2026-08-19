# Role of XIAP in Colorectal Cancer (CRC)

Schiffmann et al. (submitted) compared XIAP knock-in (KI), XIAP knock-out (KO) and wildtype mice, in order to experimentally examine the role of XIAP in CRC pathogenesis. This repository contains a reproducible workflow for the analysis of RNA-Seq data from this experiment. 

## Workflow structure

- Raw `RNA seq reads` (link to be provided) were processed with the [nf-core rnaseq pipeline v3.6](https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md)
- `quant.sf files` are produced by the [salmon](https://combine-lab.github.io/salmon/) software as part of the nf-core pipeline. They contain probabilistic `estimates of individual transcript expression` 
- A local [targets](https://books.ropensci.org/targets/) pipeline computes `Differential Gene Expression (DGE)` and `Gene Set Enrichment Analysis (GSEA)` results from the quant.sf files 

## Software structure

Two pieces of local software were used in this project:
- The `targets pipeline`
- The R script `reproduce_figures.R`: It draws figures based on the DGE and GSEA objects

Both scripts come with a `renv.lock` file describing the expected environment of execution.


## Repository structure

- `_targets.R` — R script defining a [targets](https://github.com/wlandau/targets-four-minutes) pipeline for the downstream analysis of the quant.sf files
- `_targets/` — The object store of the targets pipeline. All intermediate results are stored permanently here and can be retrieved via `tar_read(object)`
- `out/` — Tabular data written by the targets pipeline
- `report/` — Analysis reports written by the targets pipeline
- `scripts/` — A general folder for R scripts. Here, it contains only `reproduce_figures.R`, which reproduces the GSEA plot, the gene set volcano plot, and the gene volcano plot, using the objects in the _targets store.
- `R/` — A general folder for internal functions
- `logs/` — Log files from pipelines


## Environment
This project uses **renv** for reproducibility.
Run `renv::restore()` to install the correct package versions.

## Version control
Git and GitHub are configured automatically by this setup script.
