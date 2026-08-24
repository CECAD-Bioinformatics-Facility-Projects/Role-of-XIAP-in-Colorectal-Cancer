library(targets)
source("R/functions.R")

source("R/functional_enrichment.R")
options(tidyverse.quiet = TRUE)

if(!dir.exists("out")) dir.create("out")

Sys.setenv(ANNOTATION_HUB_CACHE = "/home/rstudio/project/.cache")

tar_option_set(packages = c(
	"dplyr",
	"tibble",
	"readr",
	"purrr",
	"tidyr",
	"qs",

	"AnnotationHub",
	"ensembldb",
	"AnnotationDbi",
	"GenomicFeatures",
	"tximport",
	"DESeq2",
	"ashr",
	"apeglm",
	"DOSE",
	"clusterProfiler",

	"flexdashboard",
	"ggplot2",
	"ggnewscale",
	"plotly",
	"crosstalk",
	"DT",

	# remotes::install_github("RichardJActon/colourScaleR")
	"colourScaleR"
),
debug="dge_report"
)


list(
  tar_target( ## <<<------ SET ENSEMBL VERSION TO BE USED
    use_ensembl,
    _____ENSVERSION_____
  ),

	# Differential Expression Analysis ----
	## Design Info ----
	tar_target(
		design_file,
		"design.csv",
		format = "file"
	),
	tar_target(
		design,
		readr::read_csv(design_file, col_types = "cdccc")
	),
	tar_target(padj_display_cutoff, 0.05),
	tar_target(abslfc_display_cutoff, 0.25),
	tar_target(mincov_display_cutoff, 10),
	tar_target(mincov_Nsamp_display_cutoff, 6),
	## Annotation Info ----
	tar_target(
		organism_name, "Mus musculus"
	),
	tar_target(
		txdb,
		get_enseml_txdb_annohub(organism_name, use_ensembl)
	),
	tar_target(
		tx2gene,
		get_tx2gene(txdb(),
                prefix="") ## UG, for Kashkar data
	),
	## Data Read-in ----
	tar_target(
	  quant_file_paths,
	  fs::dir_ls(
	    "results_RNA-seq/star_salmon/",
	    recurse = TRUE, regexp = "quant.sf"
	  ),
	  format = "file"
	),
	tar_target(
	  files_tab,
	  nf_salmon_quant_files_table(
	    quant_file_paths,
	    c("condition","repeat"), sep_regex = ".*\\/(\\w+)_(\\w+)\\/"
	  )
	),
	tar_target(formula, "~condition"),
	tar_target(
		quant_files,
		nf_salmon_quant_files(files_tab)
	),
	tar_target(
		txi,
		named_tximport(files_tab, type = "salmon", tx2gene = tx2gene)
	),
	tar_target(
		raw_counts,
		tximport::tximport(
			quant_files, type = "salmon", tx2gene = tx2gene,
			countsFromAbundance = "no"
		),
		format = "qs"
	),
	tar_target(
		vst_counts, DESeq2::vst(round(raw_counts$counts)), format = "qs"
	),
	tar_target(
		counts,
		raw_counts$counts %>%
			as.data.frame() %>%
			tibble::rownames_to_column("gene_id") %>%
			dplyr::left_join(
				tx2gene %>%
					tibble::as_tibble() %>%
					distinct(gene_id = GENEID, SYMBOL),
				by = "gene_id"
			) %>%
			dplyr::select(gene_id, SYMBOL, everything()),
		format = "qs"
	),
	tar_target(
		vst_counts_anno,
		vst_counts %>%
			as.data.frame() %>%
			tibble::rownames_to_column("gene_id") %>%
			dplyr::left_join(
				tx2gene %>%
					tibble::as_tibble() %>%
					distinct(gene_id = GENEID, SYMBOL),
				by = "gene_id"
			) %>%
			dplyr::select(gene_id, SYMBOL, everything()),
		format = "qs"
	),
	tar_target(
		count_matrix,
		(function(path){
			fs::dir_create(fs::path_dir(path))
			vroom::vroom_write(counts, path, delim = "\t")
			return(path)
		##})("report/count_matrix_raw.tsv"),
		})("count_matrix_raw______ENSVERSION_____.tsv"), ## UG Aug 24, 2026

		format = "file"
	),
	tar_target(
		counts_vst_file,
		(function(path){
			fs::dir_create(fs::path_dir(path))
			vroom::vroom_write(vst_counts_anno, path, delim = "\t")
			return(path)
		##})("report/count_matrix_vst.tsv"),
		})("count_matrix_vst______ENSVERSION_____.tsv"), ## UG Aug 24, 2026
		format = "file"
	),
	tar_target(
		sample_table,
		make_DESeq_sample_table(
			files_tab, as.formula(formula),
			pair = c("WT", "KO", "KI"), name = condition ## Kashkar
		)
	),
	tar_target(
		pairs,
		sample_table %>%
			dplyr::pull(condition) %>%
			levels() %>%
			combn(2, simplify = FALSE)
	),
	tar_target(
		sample_table_sub,
		sample_table %>% dplyr::filter(condition %in% pairs[[1]]),
		pattern = map(pairs)
	),

	tar_target(
		dds_txi,
		DESeq2::DESeqDataSetFromTximport(
			txi_sub(txi, col = rownames(sample_table_sub)),
			sample_table_sub, ~condition
		),
		pattern = map(sample_table_sub),
		format = "qs"
	),
	## DESseq2 differential expression ----
	tar_target(
		dds, DESeq2::DESeq(dds_txi),
		pattern = map(dds_txi),
		format = "qs"
	),
	tar_target(
		res, DESeq2::results(dds),
		pattern = map(dds),
		format = "qs"
	),
	# pairwise differential gene expression between conditions was computed
	# with DESeq2, log fold change values have been subject to shrinkage with
	# the 'ashr' method.
	tar_target(
		resLFC,
		DESeq2::lfcShrink(
			dds, coef = paste(
				"condition", pairs[[1]][2],"vs", pairs[[1]][1], sep = "_"
			), type = "ashr" #"apeglm"
		),
		pattern = map(dds, pairs),
		format = "qs"
	),
	## Result formatting -----
	tar_target(
		results_annotated,
		add_gene_annotation(resLFC, txdb()),
		pattern = map(resLFC),
		format = "qs"
	),
	### complete results tables TSVs----
	tar_target(
		results_annotated_files,
		(function(path){
			fs::dir_create(fs::path_dir(path))
			vroom::vroom_write(results_annotated, path, delim = "\t")
			return(path)
		})(paste0(
			"out/dge_results/dge_results-",
			paste0(rev(pairs[[1]]), collapse = "_vs_"),
			".tsv"
		)),
		pattern = map(results_annotated, pairs),
		format = "file"
	),
	### Min cover threshold results tables excel ----
	tar_target(
		results_annotated_min_cov,
		# coverage of at least 10 reads in at least 6 samples
		min_coverage(
			results_annotated, dds,
			min_cov = mincov_display_cutoff,
			min_samples = mincov_Nsamp_display_cutoff
		), #
		pattern = map(results_annotated, dds),
		format = "qs"
	),
	tar_target(
		results_annotated_min_cov_xl,
		results_xlsx(
			results_annotated_min_cov,
			paste0(rev(pairs[[1]]), collapse = "_vs_"),
			paste0(
				"out/dge_results/dge_results_min-cov_",
				paste0(rev(pairs[[1]]), collapse = "_vs_"),
				".xlsx"
			)
		),
		pattern = map(results_annotated_min_cov, pairs),
		format = "file"
	),
	## Combined DGE and counts spreadsheet ----
	tar_target(
		results_annotated_min_cov_grp,
		results_annotated_min_cov %>%
			dplyr::mutate(
				comparison = paste0(rev(pairs[[1]]), collapse = "_vs_")
			),
		pattern = map(results_annotated_min_cov, pairs)
	),
	tar_target(
		results_annotated_min_cov_wide,
		get_res_annot_wide(results_annotated_min_cov_grp)
	),
	tar_target(
		results_annotated_min_cov_wide_counts,
		dplyr::left_join(
			results_annotated_min_cov_wide,
			counts,
			by = "gene_id"
		)
	),
	tar_target(
		results_annotated_min_cov_wide_counts_xl,
		res_annot_wide_counts_xl(
			results_annotated_min_cov_wide_counts,
			"DGE comparisons with counts",
			paste0("out/dge_results/results_annotated_min_cov_wide_counts",
			       ".xlsx")
		)
	),
	tar_target(
		results_annotated_min_cov_wide_counts_vst,
		dplyr::left_join(
			results_annotated_min_cov_wide,
			vst_counts_anno,
			by = "gene_id"
		)
	),
	tar_target(
		results_annotated_min_cov_wide_counts_vst_xl,
		res_annot_wide_counts_xl(
			results_annotated_min_cov_wide_counts_vst,
			"DGE comparisons with counts",
			paste0("out/dge_results/results_annotated_min_cov_wide_counts_vst",
			".xlsx")
		)
	),
	# Functional Enrichment analysis ----
	tar_target(
		orgDb, get_orgDb(organism_name)
	),
	tar_target(
		ontology, c("BP", "CC", "MF")
	),

	# GSEA ----
	tar_target(
		GO_GSEA,
		dds_res_GSEA(
			results_annotated_min_cov,
			orgDb(),
			ontology = ontology,
			seed = TRUE
		),
		pattern = cross(results_annotated_min_cov, ontology),
		format = "qs"
	),
	tar_target(
		GO_GSEA_names,
		c(
			"Gene Set Enrichment, GO: Biological Processes",
			"Gene Set Enrichment, GO: Cellular Component",
			"Gene Set Enrichment, GO: Molecular Function"
		)
	),
	tar_target(
		GO_GSEA_names_exp,
		rep(GO_GSEA_names, length(pairs))
	),
	tar_target(
		GO_GSEA_names_pairs,
		get_pair_names(pairs, GO_GSEA_names)
	),
	tar_target(
		GO_GSEA_csvs,
		write_FE(
			GO_GSEA, GO_GSEA_names_pairs, "GSEA",
			"out/GSEA"
		),
		pattern = map(GO_GSEA_names_pairs, GO_GSEA),
		format = "file"
	),
	# tar_target(
	# 	GO_GSEA_res_tables_prep,
	# 	GSEA_res_table_prep(GO_GSEA),
	# 	pattern = map(GO_GSEA),
	# 	format = "qs"
	# ),
	# tar_target(
	# 	GO_GSEA_res_tables_formatted,
	# 	GO_GSEA_RT(GO_GSEA_res_tables_prep),
	# 	pattern = map(GO_GSEA_res_tables_prep),
	# 	format = "qs"
	# ),
	tar_target(
		GO_GSEA_emapplot,
		emapplot(GO_GSEA[[1]], GO_GSEA_names_exp),
		pattern = map(GO_GSEA_names_exp, GO_GSEA),
		format = "qs"
	),
	tar_target(
		GO_GSEA_gseaplot2,
		gseaplot2(GO_GSEA[[1]], GO_GSEA_names_exp, geneSetID = 1:10),
		pattern = map(GO_GSEA_names_exp, GO_GSEA),
		format = "qs"
	),
	tar_target(
		GO_GSEA_dotplot,
		dotplot(GO_GSEA[[1]], GO_GSEA_names_exp),
		pattern = map(GO_GSEA_names_exp, GO_GSEA),
		format = "qs"
	),

	# ORA ----

	## ORA up ----
	tar_target(
		GO_ORA_UP,
		dds_res_ORA(
			results_annotated_min_cov %>% filter(log2FoldChange > 0),
			orgDb(),
			ontology = ontology
		),
		pattern = cross(results_annotated_min_cov, ontology),
		format = "qs"
	),
	tar_target(
		GO_ORA_UP_names,
		c(
			"Upregulated Over-representation, GO: Biological Processes",
			"Upregulated Over-representation, GO: Cellular Component",
			"Upregulated Over-representation, GO: Molecular Function"
		)
	),
	tar_target(
		GO_ORA_UP_names_exp,
		rep(GO_ORA_UP_names, length(pairs))
	),
	tar_target(
		GO_ORA_UP_names_pairs,
		get_pair_names(pairs, GO_ORA_UP_names)
	),
	tar_target(
		GO_ORA_UP_csvs,
		write_FE(
			GO_ORA_UP, GO_ORA_UP_names_pairs, "ORA",
			"out/ORA_up"
		),
		pattern = map(GO_ORA_UP_names_pairs, GO_ORA_UP),
		format = "file"
	),

	# tar_target(
	# 	GO_ORA_UP_res_tables_prep,
	# 	ORA_res_table_prep(GO_ORA_UP),
	# 	pattern = map(GO_ORA_UP),
	# 	format = "qs"
	# ),
	# tar_target(
	# 	GO_ORA_UP_res_tables_formatted,
	# 	GO_ORA_RT(GO_ORA_UP_res_tables_prep),
	# 	pattern = map(GO_ORA_UP_res_tables_prep),
	# 	format = "qs"
	# ),
	tar_target(
		GO_ORA_UP_emapplot,
		emapplot(GO_ORA_UP[[1]], GO_ORA_UP_names_exp),
		pattern = map(GO_ORA_UP_names_exp, GO_ORA_UP),
		format = "qs"
	),
	tar_target(
		GO_ORA_UP_dotplot,
		dotplot(GO_ORA_UP[[1]], GO_ORA_UP_names_exp),
		pattern = map(GO_ORA_UP_names_exp, GO_ORA_UP),
		format = "qs"
	),

	## ORA down ----

	tar_target(
		GO_ORA_DOWN,
		dds_res_ORA(
			results_annotated_min_cov %>% filter(log2FoldChange < 0),
			orgDb(),
			ontology = ontology
		),
		pattern = cross(results_annotated_min_cov, ontology),
		format = "qs"
	),
	tar_target(
		GO_ORA_DOWN_names,
		c(
			"Downregulated Over-representation, GO: Biological Processes",
			"Downregulated Over-representation, GO: Cellular Component",
			"Downregulated Over-representation, GO: Molecular Function"
		)
	),
	tar_target(
		GO_ORA_DOWN_names_exp,
		rep(GO_ORA_DOWN_names, length(pairs))
	),
	tar_target(
		GO_ORA_DOWN_names_pairs,
		get_pair_names(pairs, GO_ORA_DOWN_names)
	),
	tar_target(
		GO_ORA_DOWN_csvs,
		write_FE(
			GO_ORA_DOWN, GO_ORA_DOWN_names_pairs, "ORA",
			"out/ORA_down"
		),
		pattern = map(GO_ORA_DOWN_names_pairs, GO_ORA_DOWN),
		format = "file"
	),

	# tar_target(
	# 	GO_ORA_DOWN_res_tables_prep,
	# 	ORA_res_table_prep(GO_ORA_DOWN),
	# 	pattern = map(GO_ORA_DOWN),
	# 	format = "qs"
	# ),
	# tar_target(
	# 	GO_ORA_DOWN_res_tables_formatted,
	# 	GO_ORA_RT(GO_ORA_DOWN_res_tables_prep),
	# 	pattern = map(GO_ORA_DOWN_res_tables_prep),
	# 	format = "qs"
	# ),
	tar_target(
		GO_ORA_DOWN_emapplot,
		emapplot(GO_ORA_DOWN[[1]], GO_ORA_DOWN_names_exp),
		pattern = map(GO_ORA_DOWN_names_exp, GO_ORA_DOWN),
		format = "qs"
	),
	tar_target(
		GO_ORA_DOWN_dotplot,
		dotplot(GO_ORA_DOWN[[1]], GO_ORA_DOWN_names_exp),
		pattern = map(GO_ORA_DOWN_names_exp, GO_ORA_DOWN),
		format = "qs"
	) ######,

	# Report formatting ----

# 	tar_target(
# 		results_formatted,
# 		format_results(results_annotated, organism_name),
# 		pattern = map(results_annotated),
# 		format = "qs"
# 	),
#
# 	tar_target(
# 		results_formatted_filt,
# 		results_formatted %>%
# 			dplyr::filter(
# 				padj < padj_display_cutoff,
# 				abs(log2FoldChange) > abslfc_display_cutoff
# 			),
# 		pattern = map(results_formatted),
# 		format = "qs"
# 	),
# 	tar_target(
# 		good_results, ## using this as index for slice() does NOT work
# 		##which(nrow(results_formatted_filt) > 1),
# 		nrow(results_formatted_filt) > 1,
# 		pattern = map(results_formatted_filt),
# 		format = "qs"
# 	),
#
# 	tar_target(
# 		shared_data,
# 		SharedData$new(results_formatted_filt),
# 		pattern = map(results_formatted_filt),
# 		format = "qs"
# 	),
# 	tar_target(
# 		dge_report_tibble,
#
# 		tibble::tibble(
# 			sample_table_sub = list(sample_table_sub),
# 			shared_data = list(shared_data),
# 			name = paste0(rev(pairs[[1]]), collapse = "_vs_"),
# 			full_dge_file = results_annotated_files,
# 			min_cov_dge_file = results_annotated_min_cov_xl,
# 			output_file = paste0(
# 				here::here(), "/report/dge_report-", name, ".html"
# 			),
# 			use=good_results ##UG
# 		),
# 		pattern = map(
# 			shared_data, pairs, sample_table_sub,
# 			results_annotated_files, results_annotated_min_cov_xl, ##?nonsense
# 			good_results ## ?nonsense
# 		),
# 		format = "qs"
# 	),
#
# 	 tarchetypes::tar_render_rep(
# 	 	dge_report, "dge_report.Rmd",
# 	 	params = dge_report_tibble %>% dplyr::filter(use) %>% dplyr::ungroup()
# 	 	#batches = 4
# 	 ) ,
# 	# Reports ----
# 	 tar_target(
# 	 	fe_report_tibble, ## suffix _full by UG
# 	 	tibble::tibble(
# 	 		ont = rep(ontology, length(pairs)),
#
# 	 		GO_GSEA_csvs = GO_GSEA_csvs,
# 	 		GO_GSEA_res_tables_formatted = GO_GSEA_res_tables_formatted,
# 	 		GO_GSEA_emapplot = GO_GSEA_emapplot,
# 	 		GO_GSEA_gseaplot2 = GO_GSEA_gseaplot2,
# 	 		GO_GSEA_dotplot = GO_GSEA_dotplot,
#
# 	 		GO_ORA_UP_csvs = GO_ORA_UP_csvs,
# 	 		GO_ORA_UP_res_tables_formatted = GO_ORA_UP_res_tables_formatted,
# 	 		GO_ORA_UP_emapplot = GO_ORA_UP_emapplot,
# 	 		GO_ORA_UP_dotplot = GO_ORA_UP_dotplot,
#
# 	 		GO_ORA_DOWN_csvs = GO_ORA_DOWN_csvs,
# 	 		GO_ORA_DOWN_res_tables_formatted = GO_ORA_DOWN_res_tables_formatted,
# 	 		GO_ORA_DOWN_emapplot = GO_ORA_DOWN_emapplot,
# 	 		GO_ORA_DOWN_dotplot = GO_ORA_DOWN_dotplot,
#
# 	 		name = purrr::map_chr(
# 	 			rep(pairs, each = length(ontology)),
# 	 			~paste0(rev(.x), collapse = "_vs_")
# 	 		),
# 	 		output_file = purrr::map_chr(
# 	 			rep(pairs, each = length(ontology)),
# 	 			~paste0(
# 	 				here::here(), "/report/fe_report-",
# 	 				paste0(rev(.x), collapse = "_vs_"), ".html"
# 	 			)
# 	 		)
# 	 	) %>%
# 	 		dplyr::select(-ont) %>%
# 	 		dplyr::group_by(name, output_file) %>%
# 	 		tidyr::nest(
# 	 			GO_ORA_UP = tidyselect::starts_with("GO_ORA_UP"),
# 	 			GO_ORA_DOWN = tidyselect::starts_with("GO_ORA_DOWN"),
# 	 			GO_GSEA = tidyselect::starts_with("GO_GSEA")
# 			),
# 		format = "qs"
# 	),
# 	tarchetypes::tar_render_rep(
#   	    fe_report, "fe_report.Rmd",
# 		params = fe_report_tibble#,
# 	    #iteration = "list"
#         #batches = 4
# 	),
#	tar_render(methods, "methods.Rmd")
)
