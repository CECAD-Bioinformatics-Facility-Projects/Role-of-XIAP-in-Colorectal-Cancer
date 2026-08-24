# Role of XIAP in Colorectal Cancer (CRC)

Schiffmann et al. (submitted) compared XIAP knock-in (KI), XIAP knock-out (KO) and wildtype mice, in order to experimentally examine the role of XIAP in CRC pathogenesis. This repository contains a reproducible workflow for the analysis of RNA-Seq data from this experiment. 

## Workflow structure

- Raw `RNA seq reads` (link to be provided) were processed with the [nf-core rnaseq pipeline v3.6](https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md)
- `quant.sf files` are produced by the [salmon](https://combine-lab.github.io/salmon/) software as part of the nf-core pipeline. They contain probabilistic `estimates of individual transcript expression` 
- A local [targets](https://books.ropensci.org/targets/) pipeline computes `Differential Gene Expression (DGE)`, `Gene Set Enrichment Analysis (GSEA)`, and other high-level result objects from the quant.sf files 
- `Downstream scripts post-process and visualize` these targets objects 


## Repository content

This repository contains the R source code for reproducing Figures 3A, 3B, and S3 of the submitted manuscript.

### R environments

The code has to deal with a partial R environment incompatibility between an initial analysis in 2022 and recent code for the figures of the submitted paper. 

The initial analysis runs a [targets pipeline](https://books.ropensci.org/targets/), which produces input objects of the figures code. Both parts are compatible with `R-4.1.2` and should be run under this R version. However the minimal version requirements of some functions have changed between the initial and the recent analysis. Therefore it is not possible to unite both in a single environment without changing the original environment of the targets pipeline, which could change the input to the figures.

Here we use a manual workaround. We provide two separate renv.lock files, describing their respective environments: `renv.lock_FOR_TARGETS` for the 2022 analysis and `renv.lock_FOR_FIGURES` for the figure code. 

**This project comes with the `TARGETS` environment pre-installed.** To switch to the `FIGURES` 
environment, enter in the rstudio terminal (or on any linux-type command line):

```
cp -a renv.lock_FOR_FIGURES renv.lock
```

then enter in the R console:

```{r switch_environment}
renv::restore(clean=TRUE)
```
Answer "y" to the package update prompt. (You can ignore warnings about the renv version.) Finally `restart R`.

To switch back to the `TARGETS` environment, copy `renv.lock_FOR_TARGETS` to `renv.lock` accordingly.

### The Targets Pipeline

There is an additional complexity here:

> In the initial analysis of this data, the pipeline had accidentally been run against the `[Ensembl](https://www.ensembl.org/) version 103` mouse transcript annotation, while the reads had actually been mapped against the version 105 genome build.  The `targets.R` script reproduces this initial run. 

This has but a minor impact on the figures, but we do offer the possibility here to produce both the Ensembl 103 and 105 versions of the figures. 

With the `TARGETS` environment activated, run

```{r run_targets103}
library(targets)
set.seed(6733)

setEnsVersion(infile =  "scripts/_targets_template.R",
              outfile = "scripts/_targets103.R",
              version = "103")
tar_make(script="scripts/_targets103.R", 
         store="my_targets103")

```

to produce the Ensembl version 103 pipeline results. 

The code creates a script `_targets103.R` and runs it, with the output sent to store "my_targets103". The seed assures that the Gene Set Enrichment Analysis (GSEA) results are stable between runs, in spite of the  built-in random component of the algorithm. 


To re-run the pipeline with `Ensembl version 105`, run

```{r run_targets105}
library(targets)
source("R/functions.R")

set.seed(6733)
setEnsVersion(infile =  "scripts/_targets_template.R",
              outfile = "scripts/_targets105.R",
              version = "105")
tar_make(script="scripts/_targets105.R", 
         store="my_targets105")

```


### The Figures Code

First switch to the FIGURES environment as described above.

To create the figures with Ensembl version 103, run

```{r figures103}
library(targets)
source("R/functions.R")

set.seed(6733)

setTargetsStore(infile =  "scripts/reproduce_figures_template.R",
                outfile = "scripts/reproduce_figures103.R",
                use_store = "my_targets103")
source("scripts/reproduce_figures103.R")

```
This code first writes and then executes a script `reproduce_figures103.R` in the `scripts` folder. This script reads its input from the targets store `my_targets103` (by changing the `use_store` argument, it can be set to use any other targets store). 

The script produces Figures 3A, 3B, and S3A of the submitted paper. The script returns them invisibly, as variables in the environment. They are printed in the Plots pane by typing their name:

Entering in the R console

- `gene_set_volcano` prints **Figure 3A** (the Gene Ontology (GO) term volcano plot)

- [`gene_volcano` is an R list with two entries]

    - `gene_volcano$min_cov` prints **Figure 3B** (the gene volcano plot)                  
      (considering only genes with at least 10 counts in at least 6 samples)
      
    - `gene_volcano$all` prints a gene volcano plot of all genes
    
- `GSEA_plot` prints **Figure S3A** (the GSEA plot) 


All figures refer to the comparison of the KI genotype versus WT.

Save the variables before you run the script again, to prevent over-writing:

```{r save_figures103}
GSEA_plot103 <- GSEA_plot
saveRDS(GSEA_plot103, file="GSEA_plot103.rds")

gene_set_volcano103 <- gene_set_volcano
saveRDS(gene_set_volcano103, file="gene_set_volcano103.rds")

gene_volcano103 <- gene_volcano
saveRDS(gene_volcano103, file="gene_volcano103.rds")
```

To re-create e.g. variable GSEA_plot103 in a later session, use
```{r extract_example}
GSEA_plot103 <- readRDS("GSEA_plot103.rds")
```

To produce the same figures with output of the pipeline using v105, run

```{r figures105}
library(targets)
source("R/functions.R")

set.seed(6733)

setTargetsStore(infile =  "scripts/reproduce_figures_template.R",
                outfile = "scripts/reproduce_figures105.R",
                use_store = "my_targets105")
source("scripts/reproduce_figures105.R")

```
This will produces the same variables as above, but now based on the Ensembl v105 targets. Save them from over-writing:

```{r save_figures105}
GSEA_plot105 <- GSEA_plot
gene_set_volcano105 <- gene_set_volcano
gene_volcano105 <- gene_volcano
```

Save the variables before you run the script again, to prevent over-writing:

```{r save_figures103}
GSEA_plot105 <- GSEA_plot
saveRDS(GSEA_plot105, file="GSEA_plot105.rds")

gene_set_volcano105 <- gene_set_volcano
saveRDS(gene_set_volcano105, file="gene_set_volcano105.rds")

gene_volcano105 <- gene_volcano
saveRDS(gene_volcano105, file="gene_volcano105.rds")
```


The differences between the results from the two Ensembl versions are minor.

