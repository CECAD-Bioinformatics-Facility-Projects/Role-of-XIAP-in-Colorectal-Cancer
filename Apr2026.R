library(targets)
library(clusterProfiler)
library(EnhancedVolcano)


source("gseaplot2_local.R")
library(ggplot2)

GO_GSEA <- setNames(tar_read(GO_GSEA),
					tar_read(GO_GSEA_names_pairs)
					)
N <- 10
gg <-
	GO_GSEA[["Gene_Set_Enrichment-GO_Biological_Processes-KI_vs_WT"]]

## ---------------- For the GSEA plot: ---------------------------------------
top_terms <- gg$Description[1:7] ## these terms share the same maximal p.adjust

## Request by Hamid Kashkar:
requested_terms <- c(setdiff(top_terms, "ribonucleoprotein complex biogenesis"),
                     gg$Description[gg$Description %in% c("cell killing", "innate immune response in mucosa")
                                   ]
                    )
# > requested_terms
# [1] "response to molecule of bacterial origin"
# [2] "humoral immune response"
# [3] "antimicrobial humoral response"
# [4] "response to lipopolysaccharide"
# [5] "antimicrobial humoral immune response mediated by antimicrobial peptide"
# [6] "cellular response to biotic stimulus"
# [7] "innate immune response in mucosa"
# [8] "cell killing"
## Have the terms ordered such that the requested terms appear on top,
## in the order of their NES [? and p.adjust] values
use_levels <- c(gg$Description[ i <- which(gg$Description %in% requested_terms)],
				gg$Description[-i]
)
use_gg <- gg
use_gg@result$Description <- factor(use_gg@result$Description,
									levels=use_levels)

use_terms <- use_gg@result$ID[use_gg@result$Description %in% requested_terms]

## Run my modification of clusterProfiler::gseaplot2(), requesting to return
## the 3 panel sub-plots individually as a list, instead of one merged plot:
p <- gseaplot2_local( ## see gseaplot_local.R
	use_gg,
	geneSetID=use_terms,
	subplots=1:3,
	raw_plotlist = TRUE,
	linewidth = 0.5
)

## Now we can augment the sub-plots, before merging them in the end.
## Remove the legend in the line plot, and set the plot borders to a weaker
## color (grey instead of black), to improve readability of inline headers:
border_theme <-
	theme(legend.position = "none",
		  panel.border = element_rect(fill=NA,
		  		  	      linetype=1,
		  		              size=0.5,#1.2,
		  			      color="grey"),
		  axis.line = element_line(colour = 'grey', size = 0.5)
	)

## This function from https://github.com/tidyverse/ggplot2/issues/4297
## allows to target the extreme ends of the x and y plot axes for placing
## a text label (without knowing about the numerical ranges of the axes).
## x=Inf, y=Inf targets the far ends of both the x and the y axis, that is,
## the upper right corner.
place_label <- function(label, size = 5, ...) {
	##	annotate("text", label = label, x = -Inf, y = Inf, ## in original function
	annotate("text", label = label, x = Inf, y = Inf,

			 ## hjust = 0, vjust = 1, size = size, ...) ## in original function
			 hjust = 1, vjust = 1, size = size, ...)

}

## Apply the modifications ...
pp <-
	c(list(
		p[[1]] +
			labs(x=NULL, y=NULL) + border_theme +
			place_label(label="Running Enrichment Score",
						size=unit(4, "pt") )

	),
	list(
		p[[3]] +
			labs(x=NULL, y=NULL) + border_theme +
			place_label(label="Ranked Log2FoldChange Values",
						size=unit(4, "pt") )
	),
	list(
		p[[2]] + border_theme
	)
	)
## .. and merge the 3 sub-plots
ppa <- aplot::plot_list(gglist = pp, ncol = 1, heights = c(0.35,0.15,0.5))

## Print the merged plot:
ppa

## works after renv::install("thomasp85/patchwork")
## but:
# Warning message:
# Removed 145480 rows containing missing values or values outside the scale range
# (`geom_label()`).

## ---------------- Gene Set Volcano Plot: -----------------------------------
library(EnhancedVolcano)
EnhancedVolcano(use_gg |> as.data.frame(),
				x="NES",
				y="pvalue", ##"p.adjust",
				lab=as.character(use_gg$Description),
				title = NULL,
				subtitle = NULL,

				boxedLabels = FALSE,
				selectLab = requested_terms,
				drawConnectors = TRUE,
				labSize = 3, #4,
				col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
				shape=21, pCutoff=1e-5,
				widthConnectors = 0.2,
				pointSize = 5,
				ylim=c(0,12), #7.5),
				colAlpha = 1,
				xlab= "Normalized Enrichment Score",
				ylab = bquote(~-Log[10] ~ pvalue))
## export as png via plot pane, setting xdim == ydim



## Gene Volcano:
library(GO.db)
library(org.Mm.eg.db)

GOALL_0002237 <- (AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0002237"), columns = c('SYMBOL'), keytype = "GOALL"))$SYMBOL |> unique()
#length(GOALL_0002237)
#[1] 409
GO_0002237    <- (AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0002237"), columns = c('SYMBOL'), keytype = "GO"))$SYMBOL |> unique()
#length(GO_0002237)
#[1] 11

## NOTE that the per-comparison sub-object of results_annotated_min_cov are
## not correctly distinguished by their listing order in _targets/objects!
## For the full dataset, _targets/objects/results_annotated_ee3a0030 is second
## in the file list and it does reproduce the second (KI_vs_WT) volcano plot
## as expected. However,
## for the filtered dataset, the object listed second,
## _targets/objects/results_annotated_min_cov_a35cca49,
## results in a completely different plot of mostly non-significant genes.
## Instead, the object listed last,
## _targets/objects/results_annotated_min_cov_920d1df3,
## produces a Volcano plot nearly identical to the plot with the full dataset.
##
## Conversely, _targets/objects/results_annotated_min_cov_a35cca49 yields a
## plot (nearly? completely?) identical to the plot of
## _targets/objects/results_annotated_c8993794 (listed third).
##
## This establishes the pairwise equivalence of one full with one filtered
## pairwise comparison at a time, but not necessarily the identity of those
## comparisons!
##
## Known equivalences:
## results_annotated_ee3a0030 == results_annotated_min_cov_920d1df3
## results_annotated_c8993794 == results_annotated_min_cov_a35cca49
## Equivalence implied by these relationships:
## results_annotated_451f5a7d == results_annotated_min_cov_dd50b015
##
## The identity of a comparisons can be verified via the resLFC* objects,
## which are input to the annotated tables.
## These objects list the underlying comparison in their heading.
##
## In one the two comparisons checked, gene Cxcl3 is highly significantly up,
## with both the unfiltered and the filtered counts.
## If this is the case for only one of **all three* comparisons, then this
## comparison can be identified by checking the gene's expression in the
## resLFC* tables.
## This is indeed the case.

# -------------------- added START -------------------------------------------

# tar_read(results_annotated_ee3a0030) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 2.744157
# tar_read(results_annotated_c8993794) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0660566
# tar_read(results_annotated_451f5a7d) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0005228838

# tar_read(results_annotated_min_cov_dd50b015) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0005228838
# tar_read(results_annotated_min_cov_a35cca49) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0660566
# tar_read(results_annotated_min_cov_920d1df3) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 2.744157


# -------------------- added END ---------------------------------------------

## The resLFC* tables contain Ensembl IDs only (no gene names),
## so we first need to get Cxcl3's ID:
Cxcl3_ID <-
	tar_read(results_annotated_ee3a0030) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(gene_id)
Cxcl3_ID
#[1] "ENSMUSG00000029379"

## Then get its expression in the 3 resLFC* tables:
tar_read(resLFC_5986ba57) |> as.data.frame() |>
	tibble::rownames_to_column("gene_id") |> dplyr::filter(gene_id==Cxcl3_ID)
#               gene_id baseMean log2FoldChange      lfcSE      pvalue      padj
# 1 ENSMUSG00000029379 25.19367   0.0005228838 0.02861681 0.003112994 0.8191799

tar_read(resLFC_691f1279) |> as.data.frame() |>
	tibble::rownames_to_column("gene_id") |> dplyr::filter(gene_id==Cxcl3_ID)
#               gene_id baseMean log2FoldChange     lfcSE       pvalue         padj
# 1 ENSMUSG00000029379 68.77935       2.744157 0.7769326 5.956506e-09 2.305049e-05

tar_read(resLFC_fab9a06a) |> as.data.frame() |>
	tibble::rownames_to_column("gene_id") |> dplyr::filter(gene_id==Cxcl3_ID)
#               gene_id baseMean log2FoldChange     lfcSE      pvalue      padj
# 1 ENSMUSG00000029379 81.95538      0.0660566 0.1565543 0.009573204 0.1818806


## Its header of resLFC_691f1279 confirms that it refers to the correct comparison:
tar_read(resLFC_691f1279)
# log2 fold change (MMSE): condition KI vs WT
# Wald test p-value: condition KI vs WT
# DataFrame with 53284 rows and 5 columns
## [...]

## Select the coresponding "results_annotated" tables via the Cxcl3 expression:
## (see above:)
# df_min_cov <- tar_read("results_annotated_min_cov_920d1df3") |>
#   mutate(GOALL_GO.0002237 = gene_name %in% GOALL_0002237)
# Error in `mutate()`:
#   ! Problem while computing `GOALL_GO.0002237 = gene_name %in% GOALL_0002237`.
# Caused by error:
#   ! `vec_is_vector()` is defunct.
# df_all <- tar_read("results_annotated_ee3a0030") |>
# 	mutate(GOALL_GO.0002237 = gene_name %in% GOALL_0002237)
# packageVersion("vctrs")
## [1] ‘0.7.3’
## Why do I have such a new version of it here??

## Workaround:

df_min_cov <- tar_read("results_annotated_min_cov_920d1df3")
df_min_cov$GOALL_GO.0002237 <-
  sapply(df_min_cov$gene_name, function(x) x %in% GOALL_0002237)

df_all <- tar_read("results_annotated_ee3a0030")
df_all$GOALL_GO.0002237 <-
  sapply(df_all$gene_name, function(x) x %in% GOALL_0002237)


DFs <- list(min_cov = df_min_cov, all = df_all)
##DFs.rds <- saveRDS(DFs, file="DFs.rds")

EV_genes_padj <-
	lapply(DFs,
		   function(df) {

		   	use_df <- df |> arrange(GOALL_GO.0002237)
		   	use_labs <- c(
		   		"S100a8",
		   		"Tnf",
		   		"Cxcl1",
		   		"Il1a",
		   		"Cxcl3",
		   		"Cemip",
		   		"Tnfsf10"
		   	)
		   	## have to-be-colored points LAST,
		   	## such that they are not overplotted!

		   	EnhancedVolcano(use_df,

		   					title = NULL,
		   					subtitle = NULL,
		   					##caption = NULL,

		   					x="log2FoldChange",
		   					y="padj", ##"pvalue",

		   					lab=use_df$gene_name,
		   					boxedLabels = FALSE,

		   					selectLab = use_labs,

		   					colCustom =
		   						setNames(ifelse(use_df$GOALL_GO.0002237,
		   								        "red",
		   								        "darkgrey"),
		   								 ifelse(use_df$GOALL_GO.0002237,
		   										"response to molecule of bacterial origin [GO:0002237]",
		   										 "other")
		   					    ),

		   					colAlpha = 0.5, ## = INtransparent

		   					drawConnectors = TRUE,
		   					labSize = 4,
		   					##col = NULL, #"black", #c('grey30', 'forestgreen', 'royalblue', 'red2'),
		   					##shape=21,
		   					pCutoff=0.05,
		   					widthConnectors = 1,
		   					pointSize = ifelse(use_df$gene_name %in% use_labs,
		   									   2,0.5), ##0.5,
		   					##ylim=c(0,7.5),
		   					xlab= "log2FoldChange",
		   					ylab = bquote(~-Log[10] ~ "adj. pvalue"), ##pvalue),

		   					max.overlaps = Inf,
		   					xlim = c(-1.5,3),
		   					ylim = c(0,5)) ##10))
		   	#p
		   }#,simplify=FALSE
	)


