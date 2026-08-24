file_substitute <-
  function(infile, outfile,
           in_pattern = "_____ENSVERSION_____",
           out_pattern = "") {

    sub(in_pattern, out_pattern, readLines(infile)) |> writeLines(con=outfile)
  }

setEnsVersion <-
  function(infile,outfile,version) {
    file_substitute(infile=infile,outfile=outfile,
                    in_pattern="_____ENSVERSION_____",out_pattern=version)
  }
setTargetsStore <-
  function(infile,outfile,use_store) {
    file_substitute(infile=infile,outfile=outfile,
                    in_pattern="_____USESTORE_____",out_pattern=use_store)
  }



#' make_DESeq_sample_table
#'
#' @param files_tab
#' @param formula
#' @param pair
#' @param name
#'
#' @return
#' @export
#'
#' @examples
make_DESeq_sample_table <- function(
	files_tab, formula, pair = NULL, name = NULL
) {
	if(!is.null(pair)) {
		files_tab <- files_tab %>%
			dplyr::filter({{name}} %in% pair) %>%
			dplyr::mutate({{name}} := factor({{name}}, levels = pair))
	}
	res_tab <- dplyr::select(files_tab, all.vars(formula)) %>%
		dplyr::mutate(
			dplyr::across(dplyr::all_of(all.vars(formula)), as.factor)
		) %>%
		as.data.frame()
	rownames(res_tab) <- files_tab$name
	res_tab
}


#' basic_DT_format
#'
#' @param data
basic_DT_format <- function(data){
	DT::datatable(
		data,
		#height = "100%",
		width = "100%",
		filter = "top",  # allows filtering on each column
		extensions = c(
			"Buttons",  # add download buttons, etc
			"Scroller"  # for scrolling down the rows rather than pagination
		),
		style = "bootstrap",
		class = "compact",
		options = list(
			dom = "Brtip",  # specify content (search box, etc)
			scroller = TRUE,
			scrollY = "960px",
			scrollX = "960px",
			buttons = list(
				"csv",  # download as .csv
				"excel"  # download as .xlsx
			)
		)
	)
}


#' get_enseml_txdb_annohub
#'
#' @param species Scientific name of a species
#' @param ensembl_version the version of ensemble to use
#'
get_enseml_txdb_annohub <- function(species, ensembl_version) {
	txdb <- function() {
		##AH <- AnnotationHub::AnnotationHub(cache = "~/.cache/AnnotationHub")
		AH <- AnnotationHub::AnnotationHub()

		AnnotationHub::query(
			AH, pattern = c(species, "EnsDb", ensembl_version)
		)[[1]]
	}
	return(txdb)
}


#' get_tx2gene
#'
#' @param txdb
#'
get_tx2gene <- function(txdb, prefix="transcript:") {
        k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
        tx2gene <- AnnotationDbi::select(
                txdb, k, c("GENEID", "SYMBOL"), "TXNAME"
        ) %>%
                dplyr::mutate(TXNAME = paste0(prefix, TXNAME))
}


#' nf_salmon_quant_files_table
#'
#' @param files
#' @param variables
#' @param sep_regex
#'
#' @return
#' @export
#'
#' @examples
nf_salmon_quant_files_table <- function(
  files, variables, sep_regex = NULL
) {
  if(is.null(sep_regex)) {
    sep_regex <- paste0(
      paste(rep("(\\w+)", length(variables)), collapse = "_"), "[\\/\\_]"
    )
  }

  tibble::tibble(
    file_path = files
  ) %>%
    tidyr::extract(
      file_path, into = variables, regex = sep_regex, remove = FALSE
    ) %>%
    dplyr::mutate(
      file_name = fs::path_ext_remove(fs::path_file(file_path)),
      name = fs::path_file(fs::path_dir(file_path))
    )
}

#' named_tximport
#'
#' @param files_tab
#' @param type
#' @param tx2gene
#'
#' @return
#' @export
#'
#' @examples
named_tximport <- function(files_tab, type = "salmon", tx2gene) {
	txi <- tximport::tximport(
		files_tab$file_path, type = type, tx2gene = tx2gene,
		countsFromAbundance = "no"
	)
	for (i in 1:3) {
		colnames(txi[[i]]) <- files_tab$name
	}
	txi$counts <- round(txi$counts)
	mode(txi$counts) <- "integer"
	txi
}

write_csv_path <- function(data, path) {
	fs::dir_create(fs::path_dir(path), recurse = TRUE)
	vroom::vroom_write(data, path, delim = ",")
	return(path)
}

nf_salmon_quant_files <- function(files_tab) {
	files <- files_tab$file_path
	names(files) <- files_tab$name
	return(files)
}

#' add_gene_annotation
#'
#' @param resLFC,
#' @param txdb
#'
add_gene_annotation <- function(resLFC, txdb) {
	res_extra <- resLFC %>%
		as.data.frame() %>%
		tibble::rownames_to_column(var = "gene_id") %>%
		arrange(pvalue) %>%
		left_join(
			GenomicFeatures::genes(txdb) %>% as_tibble(), by = "gene_id"
		) %>%
		dplyr::mutate(entrezid = as.character(entrezid)) %>%
		dplyr::select(-gene_id_version, gene_id)
}

format_results <- function(res_extra, organism_name) {
	#res_extra_text_nona <-
	res_extra %>%
		dplyr::mutate(
			`-log10(p-value)` = -log10(pvalue),
			#abslfc = abs(log2FoldChange),
			text = paste0(
				"Gene: <b>", symbol, "</b> (",gene_id,")\n",
				"Coords: ", seqnames, ":", format(start, big.mark = ","), "-",
				format(end, big.mark = ","), "\n",
				"log2(Fold Change): <b>", sprintf("%.3f", log2FoldChange),
				"</b>\n",
				"p-value: <b>", sprintf("%.3e", pvalue), "</b> (",
				"-log10(p): ", sprintf("%.2f", `-log10(p-value)`), ")\n",
				"type: ", gene_biotype, "\n",
				"Description: ",
				stringr::str_wrap(description, width = 40), "\n"
			),
			gene_id_link = paste0(
				# "href=https://wormbase.org/species/c_elegans/gene/",
				"<a target='_blank' href=",
				"https://www.ensembl.org/", gsub(" ", "_", organism_name),
				"/Gene/Summary?db=core;g=", gene_id, ">", gene_id, "</a>"
			)
		) %>%
		dplyr::select(gene_id_link, dplyr::everything()) %>%
		tidyr::drop_na(pvalue)
}

#' plotly_volcano
#'
#' @param data
#'
plotly_volcano <- function(data) {
	data %>%
		plotly::plot_ly(
			x = ~log2FoldChange, y = ~`-log10(p-value)`,
			color = ~log2FoldChange,
			size = ~baseMean,
			colors = scico::scico(9, direction = -1, palette = "roma"), #romaO
			type = 'scatter', mode = "markers",
			hoverinfo = 'text',
			text = ~text#,
			# width = "100%",
			# height = "100%"
		) %>%
		plotly::layout(
			title = "Differentially Expressed Genes (DESeq2)",
			xaxis = list(title = "log2(Fold Change)"),
			yaxis = list(title = "-log10(p-value)")#,
		) %>%
		plotly::toWebGL()
}


lfcpal <- function(x) {
	if(min(x)==max(x)) return(NULL) ## UG

	pal <- colourScaleR::universal_colour_scaler(
		x, #log2FoldChange,
		type = "scico", palette = "romaO",
		mode = "palette", direction = -1, n_breaks = 9
	)
	DT::styleInterval(as.numeric(names(pal))[-1], pal)
}

pvalpal <- function(x) {
	if(min(x)==max(x)) return(NULL) ## UG

	# pal <- colourScaleR::universal_colour_scaler(
	# 	x, #pvalue.
	# 	type = "scico", palette = "romaO",
	# 	mode = "palette", direction = -1, n_breaks = 9
	# )
	pal <- colourScaleR::universal_colour_scaler(
		x, #pvalue,
		type = "scico", palette = "hawaii", mode = "palette",
		direction = 1, n_breaks = 9#, verbose = TRUE
	)
	# pal
	DT::styleInterval(as.numeric(names(pal))[-1], pal)
}

results_table_DT <- function(data) {
	data %>%
		DT::datatable(
			filter = "top",  # allows filtering on each column
			extensions = c(
				"Buttons",  # add download buttons, etc
				"Select",
				"Scroller"  # for scrolling down the rows rather than pagination
			),
			selection = 'none',
			rownames = FALSE,  # remove rownames
			style = "bootstrap",
			class = "compact",
			width = "100%",
			height = "100%",
			escape = FALSE,
			options = list(
				dom = "Blrtip",  # specify content (search box, etc)
				select = list(style = 'os', items = 'row'),
				deferRender = TRUE,
				scrollY = "960px",
				scrollX = "960px",
				scroller = TRUE,
				columnDefs = list(
					list(
						visible = FALSE,
						#targets = c(7,8,9,10,13,15,16,17,18)
						targets = c(1, c(7,8,9,10,13,15,16,17,18) + 1)
					)
				),
				searchCols = c(
					vector(mode = "list", length = 6),
					list(list(search = '0.000 ... 0.050')),
					vector(mode = "list", length = 12)
				),
				buttons = list(
					I("colvis"),  # turn columns on and off
					# "selectRows",
					"selectNone",
					"csv",  # download as .csv
					"excel"  # download as .xlsx
				)
			)
		) %>%
		DT::formatSignif(
			c("baseMean", "pvalue", "padj", "log2FoldChange", "lfcSE"),
			digits = 4
		) %>%
		DT::formatStyle(
			"log2FoldChange",
			backgroundColor = lfcpal(.$x$data$log2FoldChange)
		) %>%
		DT::formatStyle(
			"pvalue",
			color = "#FFFFFF",
			backgroundColor = pvalpal(.$x$data$pvalue)
		) %>%
		DT::formatStyle(
			"padj",
			color = "#FFFFFF",
			backgroundColor = pvalpal(.$x$data$padj)
		)
}


txi_sub <- function(txi, col = NULL, row = NULL) {
	if(is.null(col)) {
		col <- colnames(txi$abundance)
	}
	if(is.null(row)) {
		row <- rownames(txi$abundance)
	}
	txi$abundance <- txi$abundance[row, col]
	txi$counts <- txi$counts[row, col]
	txi$length <- txi$length[row, col]
	return(txi)
}

# Excel outputs ----

#' results_xlsx
#'
#'
#' @param results_annotated annotated deseq2 results table
#' @param name name of the sheet
#' @param file output file
results_xlsx <- function(results_annotated, name, file) {

	name <- stringr::str_trunc(name, 31)

	wb <- openxlsx::createWorkbook()
	hs <- openxlsx::createStyle(textDecoration = "bold")

	openxlsx::addWorksheet(wb, sheetName = name)
	openxlsx::writeDataTable(wb, name, results_annotated, headerStyle = hs)
	openxlsx::setColWidths(
		wb, name, cols = 1:ncol(results_annotated), widths = "auto"
	)
	openxlsx::freezePane(wb, name, firstRow = TRUE)

	fs::dir_create(fs::path_dir(file))
	openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
	file
}

min_coverage <- function(
	results_annotated, dds, min_cov = 10, min_samples = 6
) {
	results_annotated %>% dplyr::filter(
		gene_id %in% names(dds[rowSums(counts(dds) >= min_cov) >= min_samples])
	)
}

#' get_res_annot_wide
#'
#'
#' @param res_annot_grp
#'
#'
get_res_annot_wide <- function(res_annot_grp) {
	res_annot_lst <-
	res_annot_grp %>%
		dplyr::select(gene_id, log2FoldChange, pvalue, comparison) %>%
		dplyr::group_by(comparison) %>%
		dplyr::group_split()

	names(res_annot_lst) <- purrr::map_chr(res_annot_lst, ~.x$comparison[1])

	res_annot_lst <- purrr::map(names(res_annot_lst), ~{
		nm <- .x
		res_annot_lst[[nm]] %>%
			dplyr::select(-comparison) %>%
			dplyr::rename_with(~paste0(nm, "_", .x), c(log2FoldChange, pvalue))
	})
	purrr::reduce(res_annot_lst, ~dplyr::full_join(.x, .y, by = "gene_id"))
}


#' res_annot_wide_counts_xl
#'
#' @param res_annot_wide_counts
#' @param name
#' @param file
#'
res_annot_wide_counts_xl <- function(res_annot_wide_counts, name, file) {
	wb <- openxlsx::createWorkbook()
	hs <- openxlsx::createStyle(textDecoration = "bold")

	openxlsx::addWorksheet(wb, sheetName = name)
	openxlsx::writeDataTable(wb, name, res_annot_wide_counts, headerStyle = hs)
	openxlsx::setColWidths(
		wb, name, cols = 1:ncol(res_annot_wide_counts), widths = "auto"
	)
	openxlsx::freezePane(wb, name, firstRow = TRUE)

	clnms <- colnames(res_annot_wide_counts)
	lfc_cols <- which(grepl("_log2FoldChange", clnms))

	pcol <- which(grepl("_pvalue", clnms))

	symcol <- which(grepl("SYMBOL", clnms))

	for (i in lfc_cols) {
		openxlsx::conditionalFormatting(
			wb, name,
			cols = i,
			type = "colourScale",
			style = c("#5e3c99", "#f7f7f7", "#e66101"),
			rows = 2:(nrow(res_annot_wide_counts) + 1)
		)
	}

	for (i in pcol) {
		openxlsx::conditionalFormatting(
			wb, name,
			cols = i,
			type = "colourScale",
			#style = c("#d01c8b", "#f7f7f7", "#4dac26"),
			style = c("#d01c8b", "#4dac26"),
			rows = 2:(nrow(res_annot_wide_counts) + 1),
			rule = c(0, 1)
			# rule = range(res_annot_wide_counts[,pcol], na.rm = TRUE)
			# rule = quantile(
			# 	res_annot_wide_counts[
			# 		2:(nrow(res_annot_wide_counts)) + 1, pcol
			# 	],
			# 	na.rm = TRUE,
			# 	probs = c(0.01, 0.99)
			# )
			# rule = range(-log10(res_annot_wide_counts[,pcol]), na.rm = TRUE)
		)
	}

	# all_hm_vals <-
	# 	res_annot_wide_counts[
	# 		2:(nrow(res_annot_wide_counts)) + 1,
	# 		(symcol + 1) : ncol(res_annot_wide_counts)
	# 	]
	#
	# all_hm_vals <- unlist(all_hm_vals)
	#
	# all_hm_vals <- all_hm_vals[!is.na(all_hm_vals)]

	# hm_col_scale <- colourScaleR::universal_colour_scaler(
	# 	all_hm_vals,
	# 	# verbose = TRUE,
	# 	scale = "quantile", type = "brewer", palette = "RdYlBu",
	# 	n_breaks = 7,
	# 	mode = "palette"
	# )

	# breaks <- as.numeric(names(hm_col_scale))

	openxlsx::conditionalFormatting(
		wb, name,
		cols = (symcol + 1) : ncol(res_annot_wide_counts),
		type = "colourScale",
		#style = c("#5e3c99", "#f7f7f7", "#e66101"),
		style = rev(c("#FC8D59", "#FFFFBF", "#91BFDB")),
		#style = c("#EF8A62", "#F7F7F7", "#67A9CF"),
		#style = hm_col_scale,
		rows = 2:(nrow(res_annot_wide_counts) + 1),
		#rule = breaks
		rule = quantile(
			res_annot_wide_counts[
				2:(nrow(res_annot_wide_counts)) + 1,
				(symcol + 1) : ncol(res_annot_wide_counts)
			],
			na.rm = TRUE,
			probs = c(0.1, 0.5, 0.99)
		)
	)

	fs::dir_create(fs::path_dir(file))
	openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
	file
}
