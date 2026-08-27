# Role of XIAP in Colorectal Cancer (CRC)

Schiffmann et al. (submitted) compared XIAP knock-in (KI), XIAP knock-out (KO) and wildtype mice, in order to experimentally examine the role of XIAP in CRC pathogenesis. This repository contains a reproducible workflow for the analysis of RNA-Seq data from this experiment. 

## Workflow structure

- Raw `RNA seq reads` (link to be provided) were processed with the [nf-core rnaseq pipeline v3.6](https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md)
- `quant.sf files` are produced by the [salmon](https://combine-lab.github.io/salmon/) software as part of the nf-core pipeline. Here, they are provided in the `sub-folders of results_RNA-seq/star_salmon/`. They contain probabilistic `estimates of individual transcript expression`. 
- A local [targets](https://books.ropensci.org/targets/) pipeline computes `Differential Gene Expression (DGE)`, `Gene Set Enrichment Analysis (GSEA)`, and other high-level result objects from the quant.sf files. 
- `Downstream scripts post-process and visualize` these targets objects. 


## Repository content

This repository contains the R source code for reproducing Figures 3A, 3B, and S3 of the submitted manuscript.

### R environments

The code has to deal with a partial R environment incompatibility between an initial analysis in 2022 and recent code for the figures of the submitted paper. 

The initial analysis runs a [targets pipeline](https://books.ropensci.org/targets/), which produces the input objects for the figures code. Both parts are compatible with `R-4.1.2` and should be run under this R version. However the minimal version requirements of some functions have changed between the initial and the recent analysis. Therefore it is not possible to unite both in a single environment without changing the original environment of the targets pipeline, which could change the input to the figures.

Here we use a manual workaround. We provide two separate renv.lock files, describing their respective environments: `renv.lock_FOR_TARGETS` for the 2022 analysis and `renv.lock_FOR_FIGURES` for the figure code. 

### The Targets Pipeline

To activate the `TARGETS` environment, enter in the R console:
```{r restore_targets}
renv::restore(lockfile = "renv.lock_FOR_TARGETS", clean=TRUE)
```
Answer "y" to the package update prompt. (You can ignore warnings about the renv version -- the next steps will take care of this.) 

Next `restart R`,  

and then enter in the console
```{r activate_targets}
renv::activate()
```

[To later switch to the `FIGURES` environment, follow these same steps, with `renv.lock_FOR_TARGETS` replaced by `renv.lock_FOR_FIGURES`.]


With the `TARGETS` environment activated, run
```{r run_targets}
library(targets)
set.seed(6733)

tar_make(script="scripts/_targets.R")
```
to produce the pipeline results in a folder named "_targets", using Ensembl version 105.  

Setting a "seed"" assures that the Gene Set Enrichment Analysis (GSEA) results produced by the pipeline are stable between runs, in spite of the  built-in random component of the algorithm. 


### The Figures Code
To switch to the `FIGURES` environment, run the same code as described above for TARGETS, but with `renv.lock_FOR_TARGETS` replaced by `renv.lock_FOR_FIGURES`.

To create the figures (using Ensembl version 105), run

```{r figures105}
library(targets)
source("R/functions.R")

set.seed(6733)

source("scripts/reproduce_figures.R")
```
This script reads its input from "_targets" and produces Figures 3A, 3B, and S3A of the submitted paper. 

The figures are returned invisibly, as variables in the environment. They are printed in the Plots pane by typing their name at the console prompt:

- `gene_set_volcano` prints **Figure 3A** (the Gene Ontology (GO) term volcano plot)

- `gene_volcano$min_cov` prints **Figure 3B** (the gene volcano plot, considering only genes with at least 10 counts in at least 6 samples)
      
- `gene_volcano$all` prints a gene volcano plot of all genes
    
- `GSEA_plot` prints **Figure S3A** (the GSEA plot) 


All figures refer to the comparison of the KI genotype versus wildtype (WT).

### Saving Variables Permanently on Disk

Environment variables are lost upon re-starting the R session or quitting R.

They can be saved on disk like so:

```{r save_figures}
saveRDS(GSEA_plot, file="GSEA_plot.rds")
saveRDS(gene_set_volcano, file="gene_set_volcano.rds")
saveRDS(gene_volcano, file="gene_volcano.rds")
```

To read a stored variable back from disk, for example `GSEA_plot`, use

```{r read_back}
GSEA_plot <- readRDS("GSEA_plot.rds")
```


