# Role of XIAP in Colorectal Cancer (CRC)

Schiffmann et al. (submitted) compared XIAP knock-in (KI), XIAP knock-out (KO) and wildtype mice, in order to experimentally examine the role of XIAP in CRC pathogenesis. This repository contains a reproducible workflow for the analysis of RNA-Seq data from this experiment. 

## Workflow structure

- Raw `RNA seq reads` (link to be provided) were processed with the [nf-core rnaseq pipeline v3.6](https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md)
- `quant.sf files` are produced by the [salmon](https://combine-lab.github.io/salmon/) software as part of the nf-core pipeline. They contain probabilistic `estimates of individual transcript expression` 
- A local [targets](https://books.ropensci.org/targets/) pipeline computes `Differential Gene Expression (DGE)`, `Gene Set Enrichment Analysis (GSEA)`, and other high-level result objects from the quant.sf files 
- `Downstream scripts post-process and visualize` these targets objects 


## Repository structure

This repository has to deal with a partial R environment incompatibility between an initial analysis in 2022 and recent code for figures of the submitted paper. 

```{r test}
a <- b
```

- `reproduce_targets/`: This rstudio project folder is meant for the documentation and possible re-execution of the 
targets pipeline. It contains the _targets.R pipeline script, local functions in R/, a renv.lock file documenting the package environment, and the metadata of the experiment in `design.csv`. See the folder specific README file for instructions on how to set up the rstudio project, download the input quant.sf files from [Figshare](), and run the pipeline. 

- `reproduce_figures/`: This separate rstudio project folder contains the code for reproducing Figures x and y of the paper, plus local functions in R/ and the package environment in renv.lock. The code depends on the objects created by the targets pipeline. We provide a tarball of the pre-computed objects on [Figshare]() (alternatively you can create the objects yourself by re-running the pipeline). See the folder specific README for details.

Both projects use R-4.1.2, the R version under which the the original targets pipeline was written. Nevertheless they must be kept separate because of version conflicts for some packages.  


## Environment
This project uses **renv** for reproducibility.
Run `renv::restore()` to install the correct package versions.

## Version control
Git and GitHub are configured automatically by this setup script.
