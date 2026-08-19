# Role of XIAP in Colorectal Cancer

## Repository structure

- `results_RNA-seq/star_salmon/` — Gene expression quantification by the [salmon]() tool, as called by the [nf-core v3.6]() RNA-seq pipeline. Folders WT, KO, and KI contain the raw quantifications (quant_genes.sf and quant.sf) for genes and transcripts. In addition we provide tabular representations of the counts and different count transformations, as tab-separated text files (.tsv) or as RDS-compressed R SummarizedExperiment objects (.rds). 
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
