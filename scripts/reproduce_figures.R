library(dplyr)
library(targets)
library(clusterProfiler)
library(EnhancedVolcano)
source("~/project/R/gseaplot2_local.R")

use_store <- "/home/rstudio/project/_targets"

GO_GSEA <- setNames(tar_read(GO_GSEA, store=use_store),
		    tar_read(GO_GSEA_names_pairs, store=use_store)
	   )
## The comparison to use:
gg <- GO_GSEA[["Gene_Set_Enrichment-GO_Biological_Processes-KI_vs_WT"]]

use_gg <- gg

## ---------------- Constructing the GSEA plot: ---------------------------------------
top_terms <- gg$Description[1:7] ## these terms share the same maximal p.adjust

## Request by Hamid Kashkar:
requested_terms <- c(setdiff(top_terms, "ribonucleoprotein complex biogenesis"), ## do not consider this top term
                     gg$Description[gg$Description %in% c("cell killing", ## add these terms, although not among the maximal terms
                                                          "innate immune response in mucosa")
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

## Define the term order such that the requested terms appear to top:
use_levels <- c(gg$Description[ i <- which(gg$Description %in% requested_terms)],
		gg$Description[-i])

use_gg@result$Description <- factor(use_gg@result$Description, levels=use_levels)

use_terms <- use_gg@result$ID[use_gg@result$Description %in% requested_terms]

## --------------- Construct and draw an adapted version of clusterProfiler's GSEA plot:

## (1) Run R/gseaplot2_local.R (a modification of clusterProfiler::gseaplot2() ),
##     which allows to return the 3 panel sub-plots individually as a list,
##     rather than one merged plot:
p <- gseaplot2_local( ## see gseaplot_local.R
	use_gg,
	geneSetID=use_terms,
	subplots=1:3,
	raw_plotlist = TRUE,
	linewidth = 0.5
)

## (2) Augment the sub-plots:
## (2a) Remove the legend in the line plot, and set the plot borders to a weaker
##      color (grey instead of black), to improve readability of inline headers:
border_theme <-
	theme(legend.position = "none",
		  panel.border = element_rect(fill=NA,
		  		  	      linetype=1,
		  		              size=0.5,#1.2,
		  			      color="grey"),
		  axis.line = element_line(colour = 'grey', size = 0.5)
	)

## (2b) This function from https://github.com/tidyverse/ggplot2/issues/4297
##      allows to target the extreme ends of the x and y plot axes for placing
##      a text label (without knowing about the numerical ranges of the axes).
##      x=Inf, y=Inf targets the far ends of both the x and the y axis, that is,
##      the upper right corner.
place_label <- function(label, size = 5, ...) {
	##	annotate("text", label = label, x = -Inf, y = Inf, ## in original function
	annotate("text", label = label, x = Inf, y = Inf,

			 ## hjust = 0, vjust = 1, size = size, ...) ## in original function
			 hjust = 1, vjust = 1, size = size, ...)

}

## (2c) Apply the modifications, and re-order the list elements (= the sub-plots):
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
## (2d) Merge the 3 sub-plots
GSEA_plot <- aplot::plot_list(gglist = pp, ncol = 1, heights = c(0.35,0.15,0.5))

## (2e) To print the merged plot type:
## GSEA_plot


## --------------- Gene Set Volcano Plot: ----------------------------------------------
gene_set_volcao <-
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



## --------------- Gene Volcano Plot: --------------------------------------------------
## We want to draw a volcano plot for the genes in the dataset which belong to GO:0002237
## (response to molecule of bacterial origin), specifically for the KI vs WT comparison.
##
## There is a slight problem here: the identity of the individual comparisons is not obvious
## from the internal identifiers of the _targets objects.
##
## I solved this brute force by the following recipe:

## (1) In the original Volcano plots, gene Cxcl3 is highly differential in KI vs WT

## (2) In the results of the current pipeline, this gene is differential only in one of the 3 comparisons:
# tar_read(results_annotated_ee3a0030) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 2.744157
# tar_read(results_annotated_c8993794) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0660566
# tar_read(results_annotated_451f5a7d) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0005228838

## (3) The identity CAN be established in the resLC objects,
##     because these objects contain DESeq2 output files, which do list the comparison.
##     However, the resLFC objects list Ensembl IDs only, not gene names.
##     Therefore:

## (3a) Get the Ensembl ID of Cxcl3:
# Cxcl3_ID <-
#   tar_read(results_annotated_ee3a0030) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(gene_id)
# Cxcl3_ID
# [1] "ENSMUSG00000029379"

## (3b) Find the resLFC object in which Cxcl3 is differential:
# tar_read(resLFC_5986ba57) |> as.data.frame() |>
#   tibble::rownames_to_column("gene_id") |> dplyr::filter(gene_id==Cxcl3_ID)
#                 gene_id baseMean log2FoldChange      lfcSE      pvalue      padj
#   1 ENSMUSG00000029379 25.19367   0.0005228838 0.02861681 0.003112994 0.8191799
#
# tar_read(resLFC_691f1279) |> as.data.frame() |>
#   tibble::rownames_to_column("gene_id") |> dplyr::filter(gene_id==Cxcl3_ID)
#                 gene_id baseMean log2FoldChange     lfcSE       pvalue         padj
#   1 ENSMUSG00000029379 68.77935       2.744157 0.7769326 5.956506e-09 2.305049e-05
#
# tar_read(resLFC_fab9a06a) |> as.data.frame() |>
#   tibble::rownames_to_column("gene_id") |> dplyr::filter(gene_id==Cxcl3_ID)
#                 gene_id baseMean log2FoldChange     lfcSE      pvalue      padj
#   1 ENSMUSG00000029379 81.95538      0.0660566 0.1565543 0.009573204 0.1818806

## Query the comparison:
# tar_read(resLFC_691f1279)
#   log2 fold change (MMSE): condition KI vs WT
#   Wald test p-value: condition KI vs WT
#   DataFrame with 53284 rows and 5 columns
#   ## [...]

## This tells us that a log2FoldChange of 2.744157 in Cxcl3 indeed identifies the KI vs WT comparison.
## Now we can conclude that the results_annotated_ee3a0030 object also comes from this comparison,
## although this is not directly documented in the object itself.

## The same logic allows to identify the min_cov table corresponding to KI vs WT:
# tar_read(results_annotated_min_cov_dd50b015) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0005228838
# tar_read(results_annotated_min_cov_a35cca49) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 0.0660566
# tar_read(results_annotated_min_cov_920d1df3) |> dplyr::filter(gene_name=="Cxcl3") |> dplyr::pull(log2FoldChange)
# [1] 2.744157


## Gene Volcano:
### BiocManager::install("org.Mm.eg.db") ## not via renv() -- OK?

library(GO.db)
library(org.Mm.eg.db)

## GO:0002237 is "response to molecule of bacterial origin"

GOALL_0002237 <- (AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0002237"),
                                        columns = c('SYMBOL'),
                                        keytype = "GOALL"))$SYMBOL |> unique()
#length(GOALL_0002237)
#[1] 409
GO_0002237    <- (AnnotationDbi::select(org.Mm.eg.db, keys=c("GO:0002237"),
                                        columns = c('SYMBOL'),
                                        keytype = "GO"))$SYMBOL |> unique()
#length(GO_0002237)
#[1] 11



df_min_cov <- tar_read("results_annotated_min_cov_920d1df3", store=use_store) |>
	transform(GOALL_GO.0002237 = gene_name %in% GOALL_0002237)

df_all <- tar_read("results_annotated_ee3a0030", store=use_store) |>
	transform(GOALL_GO.0002237 = gene_name %in% GOALL_0002237)
## "mutate" replaced by "transform",
## because the former depends on the deprecated vctrs::vec_is_vector(),
## which is defunct -- however only since Apr 17, 2026!

DFs <- list(min_cov = df_min_cov, all = df_all)
DFs.rds <- saveRDS(DFs, file="DFs.rds")

gene_volcano <-
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
		   }
	)


## Printing a list element plots the respective volcano:
##gene_volcano$min_cov
##gene_volcano$all
