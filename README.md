# Role of XIAP in Colorectal Cancer (CRC)

Schiffmann et al. (submitted) compared XIAP knock-in (KI), XIAP knock-out (KO) and wildtype mice, in order to experimentally examine the role of XIAP in CRC pathogenesis. This repository contains a reproducible workflow for the analysis of RNA-Seq data from this experiment. 

## Workflow structure

- Raw `RNA seq reads` (link to be provided) were processed with the [nf-core rnaseq pipeline v3.6](https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md)
- `quant.sf files` are produced by the [salmon](https://combine-lab.github.io/salmon/) software as part of the nf-core pipeline. They contain probabilistic `estimates of individual transcript expression` 
- A local [targets](https://books.ropensci.org/targets/) pipeline computes `Differential Gene Expression (DGE)`, `Gene Set Enrichment Analysis (GSEA)`, and other high-level result objects from the quant.sf files 
- `Downstream scripts post-process and visualize` these targets objects 


## Repository structure

This repository has to deal with a partial R environment incompatibility between an initial analysis in 2022 and recent code for figures of the submitted paper. In addition, transcripts in the initial analysis had been annotated using Ensembl v103, while the reads had actually been mapped to the v105 mouse genome build. 

### R environments

The initial analysis runs a [targets pipeline](https://books.ropensci.org/targets/), the results of which are input to the figure code. Both are compatible with `R-4.1.2` and should be run under this R version. However the minimal version requirements of some functions have changed between the initial and the recent analysis. Therefore it is not possible to unite both in a single environment without changing the original environment of the targets pipeline, which could change the input to the figures.

Here we use a manual workaround. We provide two separate renv.lock files, describing their respective environments: `renv.lock_FOR_TARGETS` for the 2022 analysis and `renv.lock_FOR_FIGURES` for the figure code. To activate the `FIGURES` environment, enter on a linux-type command line:

```
cp -a renv.lock_FOR_FIGURES renv.lock
```

then restart R and enter in the R console:

```{r switch_environment}
renv::restore(clean=TRUE)
```
and re-start R again.

To switch to the `TARGETS` environment, copy `renv.lock_FOR_TARGETS` to `renv.lock accordingly.

### The Targets Pipeline

After activating the `TARGETS` environment, do 
```{r switch_environment}
library(targets)
set.seed(6733)

tar_make(script="scripts/targets.R", store=[MY_STORE])
```

where `[MY_STORE]` is a folder name. This is where the pipeline output goes. The name can be freely choosen, but it is typically `_targets`. The seed assures the Gene Set Enrichment Analysis (GSEA) results are stable between runs, in spite of the  built-in random component of the algorithm.

Running the pipeline is only necessary if you want to verify the output. If you are only interested in using the output, we offer two pre-filled target stores: `precomputed_targets103` and `precomputed_targets105`. The numbers refer to versions of the [Ensembl](https://www.ensembl.org/) database. The initial pipeline run had accidentally used the v103 mouse transcript annotation, while the reads had actually been mapped against the v105 genome build. 

The separate folders allow to judge the input of the version difference (which is minor).  



### The Figure Code


`



## Environment
This project uses **renv** for reproducibility.
Run `renv::restore()` to install the correct package versions.

## Version control
Git and GitHub are configured automatically by this setup script.
