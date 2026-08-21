get_orgDb <- function(organism_name,
                      use_cache="/home/rstudio/project/.cache/AnnotationHub" ## UG added Aug 21, 2026
                      ) {
	orgDb <- function() {
		AH <- AnnotationHub::AnnotationHub(cache = use_cache)
		OrgDB <- AnnotationHub::query(
			AH, pattern = c(organism_name, "OrgDB")
		)[[1]]
	}
	return(orgDb)
}

#'
#'
dds_res_GSEA <- function(
	dds_res_anno, orgDb, ontology, keyType = "SYMBOL", ...
) {
	dds_res_anno_lfcsort <- dds_res_anno %>%
		dplyr::arrange(desc(log2FoldChange))
	dds_res_anno_lfcsort_gse <- dds_res_anno_lfcsort$log2FoldChange
	names(dds_res_anno_lfcsort_gse) <- dds_res_anno_lfcsort$symbol

	res <- clusterProfiler::gseGO(
		geneList = dds_res_anno_lfcsort_gse,
		ont = ontology,
		OrgDb = orgDb,
		keyType = keyType,
		...
	)
	res <- list(res)
	#
	#names(res) <- ontology
	return(res)
}

##dds_res_ORA <- function(
##	dds_res_anno, orgDb, ontology, FDR = 0.05, keyType = "SYMBOL", ...
##) {
##	DE_genes <- dds_res_anno %>%
##		dplyr::filter(padj < FDR) %>%
##		dplyr::pull(symbol)
##
##	res <- clusterProfiler::enrichGO(
##		gene = DE_genes,
##		universe = dds_res_anno$symbol,
##		ont = ontology,
##		OrgDb = orgDb,
##		keyType = keyType,
##		...
##	)
##
##	res@result <- res@result %>%
##		dplyr::filter(p.adjust < res@pvalueCutoff, qvalue < res@qvalueCutoff)
##
##	res <- list(res)
##	#names(res) <- ontology
##	return(res)
##}

emptyEnrichResult <- function(ontology,
							  organism,
							  keyType) { ## UG
	## dummy empty data.frame:
	## by first setting, then removing a single row,
	## the colnames and the column types are retained
	df <- data.frame(ID="1",
					 Description="1",
					 GeneRatio="1",
					 BgRatio="1",
					 pvalue=1.0,
					 p.adjust=1.0,
					 qvalue=1.0,
					 geneID="1",
					 Count=1)[-1,]

	new("enrichResult", result = df,
		gene=c("1")[-1], ## c() does NOT work!
		ontology=ontology,
		organism=organism,
		keytype=keyType)
}
dds_res_ORA <- function(
		dds_res_anno, orgDb, ontology, FDR = 0.05, keyType = "SYMBOL", ...
) {
	DE_genes <- dds_res_anno %>%
		dplyr::filter(padj < FDR) %>%
		dplyr::pull(symbol)

	if(length(DE_genes)==0) { ## BEGIN UG
		res <- emptyEnrichResult(ontology,
								 AnnotationDbi::species(orgDb),
								 keyType)
	} else { ## END UG
		res <- clusterProfiler::enrichGO(
			gene = DE_genes,
			universe = dds_res_anno$symbol,
			ont = ontology,
			OrgDb = orgDb,
			keyType = keyType,
			...
		)
	} ## BEGIN UG
	if(is.null(res)) {
		res <- emptyEnrichResult(ontology,
								 AnnotationDbi::species(orgDb),
								 keyType)
	} else { ## END UG
		res@result <- res@result %>%
			dplyr::filter(p.adjust < res@pvalueCutoff, qvalue < res@qvalueCutoff)
	}
	res <- list(res)
	#names(res) <- ontology
	return(res)
}

GSEA_res_table_prep <- function(GSEA_result) {
	GSEA_result <- GSEA_result[[1]]@result
	if(nrow(GSEA_result) == 0) {
		GSEA_result <- dplyr::bind_cols(
			GSEA_result, tibble::tibble(
				rank = integer(), leading_edge = character(),
				core_enrichment = character()
			)
		) %>% as.data.frame()
	}

	GSEA_result %>%
		tibble::as_tibble() %>%
		dplyr::mutate(
			core_enrichment_l = strsplit(core_enrichment, split = "/"),
			n = purrr::map_int(core_enrichment_l, length),
			percent = (n / setSize) * 100
		) %>%
		dplyr::select(
			ID, Description, n, setSize, percent, enrichmentScore, NES,
			pvalue, p.adjust, qvalues, rank, core_enrichment
		)
}

ORA_res_table_prep <- function(ORA_result) {
  df <- tibble::tibble(
    ID = character(),
    Description = character(),
    k = integer(),
    M = integer(),
    percent = double(),
    GeneRatio = double(),
    n = integer(),
    BgRatio = double(),
    N = integer(),
    pvalue = double(),
    p.adjust = double(),
    qvalue = double(),
    geneID = character()
  )

  if(is.null(ORA_result[[1]])) {
    return(df)
  }

  if(nrow(ORA_result[[1]]@result) == 0) {
    return(df)
  }

  ORA_result <- ORA_result[[1]]@result %>% as.data.frame()
  ORA_result %>%
    tibble::as_tibble() %>%
    # see:
    # - https://github.com/YuLab-SMU/DOSE/blob/aeea30e817a5b9fc3056364af7a30e21b5092841/R/enricher_internal.R
    # - http://yulab-smu.top/clusterProfiler-book/chapter2.html#over-representation-analysis
    # - https://doi.org/10.1093/bioinformatics/bth456
    tidyr::separate(
      GeneRatio,
      #into = c("# DEGs in set (k)", "# DEGs in any set (n)"),
      into = c("k", "n"),
      sep = "/", convert = TRUE
    ) %>%
    tidyr::separate(
      BgRatio,
      # into = c("# Genes in set (M)", "# Genes with a GO term (N)"),
      into = c("M", "N"),
      sep = "/", convert = TRUE
    ) %>%
    dplyr::mutate(
      GeneRatio = k/n, BgRatio = M/N, percent = (k / M) * 100
    ) %>%
    dplyr::select(
      ID, Description, k, M, percent, GeneRatio, n, BgRatio, N,
      pvalue, p.adjust, qvalue, geneID
    )
}

GSEA_res_DT_table <- function(GSEA_tab) {
	DT::datatable(
		data = GSEA_tab,
		filter = 'top',
		# fillContainer = TRUE,
		extensions = c('Buttons', 'Select', 'Scroller'),
		selection = 'none',
		rownames= FALSE,
		style = "bootstrap",
		class = "compact",
		width = "100%",
		height = "100%",
		options(
			#pageLength = 25,
			dom = 'Blrtip',#dom = 'Blfrtip',
			select = list(style = 'os', items = 'row'),
			deferRender = TRUE,
			scroller = TRUE,
			buttons = list('colvis', 'selectNone', 'csv', 'excel'), #'copy',, 'pdf', 'print'
			# lengthMenu = list(c(10, 25, 50, -1), c(10, 25, 50, "All")),
			scrollX = 960,#"960px",
			scrollY = 960,#"960px",
			columnDefs = list(
				visible = FALSE,
				targets = 12 #ncol(GSEA_tab)+1 # "core_enrichment"
			)
		)#,
		#escape = FALSE# c("ens_gene")
	) %>%
		DT::formatSignif(
			c("percent", "pvalue", "p.adjust", "qvalues", "enrichmentScore", "NES"), digits = 4 #,"qvalue" !s?
		) %>%
		DT::formatStyle(
			"pvalue",
			backgroundColor = pvalpal(GSEA_tab$pvalue)
		) %>%
		DT::formatStyle(
			"p.adjust",
			backgroundColor = pvalpal(GSEA_tab$p.adjust)
		) %>%
		DT::formatStyle(
			"qvalues",
			backgroundColor = pvalpal(GSEA_tab$qvalues)
		)
}

ORA_res_DT_table <- function(ORA_tab) {
	DT::datatable(
		data = ORA_tab,
		filter = 'top',
		# fillContainer = TRUE,
		extensions = c('Buttons', 'Select', 'Scroller'),
		selection = 'none',
		rownames= FALSE,
		style = "bootstrap",
		class = "compact",
		width = "100%",
		height = "100%",
		options(
			#pageLength = 25,
			dom = 'Blrtip',#dom = 'Blfrtip',
			select = list(style = 'os', items = 'row'),
			deferRender = TRUE,
			scroller = TRUE,
			buttons = list('colvis', 'selectNone', 'csv', 'excel'), #'copy',, 'pdf', 'print'
			# lengthMenu = list(c(10, 25, 50, -1), c(10, 25, 50, "All")),
			scrollX = 960, # "960px",
			scrollY = 960, #"960px",
			columnDefs = list(
				visible = FALSE,
				targets = 12# ncol(ORA_tab)+1 # "geneID"
			)
		)#,
		#escape = FALSE# c("ens_gene")
	) %>%
		DT::formatSignif(
			c("percent", "pvalue", "p.adjust", "qvalue", "GeneRatio", "BgRatio"), digits = 4 #,"qvalue" !s?
		) %>%
		DT::formatStyle(
			"pvalue",
			backgroundColor = pvalpal(ORA_tab$pvalue)
		) %>%
		DT::formatStyle(
			"p.adjust",
			backgroundColor = pvalpal(ORA_tab$p.adjust)
		) %>%
		DT::formatStyle(
			"qvalue",
			backgroundColor = pvalpal(ORA_tab$qvalue)
		)
}

dotplot <- function(data, title, ...) {
  if(is.null(data)) { return(list("Insufficient data for plot")) }
  if(nrow(data@result) == 0) { return(list("Insufficient data for plot")) }
  suppressMessages(plot <- list(
    enrichplot::dotplot(data, orderBy = "x", ...) +
      scico::scale_colour_scico(palette = "batlow") +
      ggplot2::labs(title = title)
  ))
  plot
}


emapplot <- function(data, title, ...) {
  if(is.null(data)) { return(list("Insufficient data for plot")) }
  if(nrow(data@result) < 2) { return(list("Insufficient data for plot")) }
  suppressMessages(plot <- list(
    enrichplot::emapplot(enrichplot::pairwise_termsim(data), ...) +
      scico::scale_colour_scico(palette = "batlow") +
      ggplot2::labs(title = title)
  ))
  plot
}

gseaplot2 <- function(data, geneSetID = 1:10, ...) {
  #geneSetID = 1:10
  if(is.null(data@result)) {
    return(list("Insufficient data for plot"))
  }

  nrres <- nrow(data@result)
  if(nrres == 0) {
    return(list("Insufficient data for plot"))
  }

  if(max(geneSetID) > nrres) {
    geneSetID <- 1:nrres
  }
  #suppressMessages(
  plot <- list(enrichplot::gseaplot2(data, geneSetID = geneSetID, ...))
  #)
  plot
}

#' gen_pathview_plot_section <- function(
#' 	title, data, type, fnm = NULL, testing = FALSE
#' ) {
#' 	if(is.null(fnm)) {
#' 		fnm <- gsub("[^\\w]+", "_", title, perl = TRUE)
#' 	}
#' 	# data_tab <- ifelse(type == "GSEA",
#' 	# 	data@result %>% GSEA_res_table_prep(),# %>% GO_GSEA_RT(),
#' 	# 	data@result %>% ORA_res_table_prep() #%>% GO_ORA_RT()
#' 	# )
#' 	text <- c(
#' 		'',
#' 		'Plots {.tabset data-height="1000"}',
#' 		'----------------------------------------',
#' 		'',
#' 		paste0('__', title,'__'),
#' 		'',
#' 		'### Data Table',
#' 		'',
#' 		'```{r}',
#' 		'downloadthis::download_file(
#' 			path = paste0("out/", ', fnm, ', ".csv"),
#' 			output_name = fnm, output_extension = ".csv",
#' 			button_label = "Download .csv", button_type = "default",
#' 			has_icon = TRUE, icon = "fa fa-save"
#' 		)
#' 		',
#' 		#, self_contained = TRUE
#' 		'```',
#' 		'',
#' 		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE}',
#' 		# 'DT::datatable(data@result)',
#' 		ifelse(type == "GSEA",
#' 			   #'data_tab %>% GO_GSEA_RT()',
#' 			   #'data_tab %>% GO_ORA_RT()'
#' 			   'data@result %>% GSEA_res_table_prep() %>% GO_GSEA_RT()',
#' 			   #"",
#' 			   'data@result %>% ORA_res_table_prep() %>% GO_ORA_RT()'
#' 		),
#' 		'```',
#' 		switch(
#' 			type,
#' 			"GSEA" = c(
#' 				'',
#' 				'### GSEA plot',
#' 				'',
#' 				'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE}',
#' 				'enrichplot::gseaplot2(data, geneSetID = 1:10)', # + ggplot2::labs(title = title)
#' 				'```'
#' 			),
#' 			"ORA" = NULL
#' 		),
#' 		'',
#' 		'### Dotplot',
#' 		'',
#' 		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE}',
#' 		'dotplot(data, title, showCategory = 20)',
#' 		'```',
#' 		'',
#' 		'### Network',
#' 		'',
#' 		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE}',
#' 		'emapplot(data, title)',
#' 		'```',
#' 		''
#' 	)
#'
#' 	if (testing == TRUE) {
#' 		return(text)
#' 	}
#' 	knitr::knit_child(text = text, envir = environment(), quiet = TRUE)
#' }

gen_pathview_plot_section <- function(
	title,# data, type, fnm = NULL,
	data,
	testing = FALSE
) {
	gsea_plot <- ifelse("gseaplot2" %in% colnames(data), "GSEA", "ORA")
	text <- c(
		# '```{r}',
		# '',
		# "getwd()",
		# '',
		# '```'
		'',
		'Plots {.tabset data-height="1000"}',
		'----------------------------------------',
		'',
		paste0('__', title,'__'),
		'',
		'### Data Table',
		'',
		'```{r}',
		'downloadthis::download_file(',
			paste0(
				'	path = "',
				data %>% dplyr::select(tidyselect::ends_with("_csvs")) %>%
					dplyr::pull(),# %>% {fs::path_join(c(here::here(), .))},
				'",'
			),
			paste0(
				'	output_name = "',
				data %>%
					dplyr::select(tidyselect::ends_with("_csvs")) %>%
					dplyr::pull() %>%
					fs::path_file() %>%
					fs::path_ext_remove(),
				'",'
			),
			'	output_extension = ".csv",',
			'	button_label = "Download .csv", button_type = "default",',
			'	has_icon = TRUE, icon = "fa fa-save"',
		')',
		#'',
		#, self_contained = TRUE
		'```',
		'',
		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
		# 'DT::datatable(data@result)',
		'data %>% dplyr::select(tidyselect::ends_with("_res_tables_formatted")) %>% dplyr::pull()',
		'```',
		switch(
			gsea_plot,
			"GSEA" = c(
				'',
				'### GSEA plot',
				'',
				'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
				'data %>% dplyr::select(tidyselect::ends_with("_gseaplot2")) %>% dplyr::pull()',
				'```'
			),
			"ORA" = NULL
		),
		'',
		'### Dotplot',
		'',
		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
		#'dotplot(data, title, showCategory = 20)',
		'data %>% dplyr::select(tidyselect::ends_with("_dotplot")) %>% dplyr::pull()',
		'```',
		'',
		'### Network',
		'',
		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
		# 'emapplot(data, title)',
		'data %>% dplyr::select(tidyselect::ends_with("_emapplot")) %>% dplyr::pull()',
		'```',
		''
	)

	if (testing == TRUE) {
		return(text)
	}
	hash <- digest::digest(data)
	knitr::knit_child(
		text = text,
		envir = environment(),
		#envir = parent.frame(),
		quiet = TRUE,
		options = list(
			fig.path = paste0(knitr::opts_chunk$get("fig.path"), hash, "-"),
			cache.path = paste0(knitr::opts_chunk$get("cache.path"), hash, "-")
			# fig.path = fs::path_abs(paste0(knitr::opts_chunk$get("fig.path"), hash, "-")),
			# cache.path = fs::path_abs(paste0(knitr::opts_chunk$get("cache.path"), hash, "-"))
		)
	)
}


#' gen_pathview_plot_section <- function(
#' 	#title, data, type, fnm = NULL,
#' 	data,
#' 	testing = FALSE
#' ) {
#' 	gsea_plot <- ifelse("gseaplot2" %in% colnames(data), "GSEA", "ORA")
#' 	text <- c(
#' 		'',
#' 		'Plots {.tabset data-height="1000"}',
#' 		'----------------------------------------',
#' 		'',
#' 		#paste0('__', title,'__'),
#' 		'',
#' 		'### Data Table',
#' 		'',
#' 		'```{r}',
#' 		'downloadthis::download_file(
#' 			path = ',
#' 		data %>%
#' 			dplyr::select(tidyselect::ends_with("_csvs")) %>%
#' 			dplyr::pull(),
#' 		',
#' 			output_name = ',
#' 		data %>%
#' 			dplyr::select(tidyselect::ends_with("_csvs")) %>%
#' 			dplyr::pull() %>%
#' 			fs::path_file() %>%
#' 			fs::path_ext_remove()
#' 		,',
#' 			output_extension = ".csv",
#' 			button_label = "Download .csv", button_type = "default",
#' 			has_icon = TRUE, icon = "fa fa-save"
#' 		)
#' 		',
#' 		#, self_contained = TRUE
#' 		'```',
#' 		'',
#' 		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
#' 		# 'DT::datatable(data@result)',
#' 		data %>%
#' 			dplyr::select(tidyselect::ends_with("_res_tables_formatted")) %>%
#' 			dplyr::pull(),
#' 		'```',
#' 		switch(
#' 			gsea_plot,
#' 			"GSEA" = c(
#' 				'',
#' 				'### GSEA plot',
#' 				'',
#' 				'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
#' 				data %>%
#' 					dplyr::select(tidyselect::ends_with("_gseaplot2")) %>%
#' 					dplyr::pull(),
#' 				'```'
#' 			),
#' 			"ORA" = NULL
#' 		),
#' 		'',
#' 		'### Dotplot',
#' 		'',
#' 		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
#' 		#'dotplot(data, title, showCategory = 20)',
#' 		data %>%
#' 			dplyr::select(tidyselect::ends_with("_dotplot")) %>%
#' 			dplyr::pull(),
#' 		'```',
#' 		'',
#' 		'### Network',
#' 		'',
#' 		'```{r, out.width = "960px", out.height = "960px", dpi = 600, message=FALSE, warning=FALSE}',
#' 		# 'emapplot(data, title)',
#' 		data %>%
#' 			dplyr::select(tidyselect::ends_with("_emapplot")) %>%
#' 			dplyr::pull(),
#' 		'```',
#' 		''
#' 	)
#'
#' 	if (testing == TRUE) {
#' 		return(text)
#' 	}
#' 	knitr::knit_child(text = text, envir = environment(), quiet = TRUE)
#' }

bar_chart <- function(
	label, width = "100%", height = "14px", fill = "#00bfc4", background = NULL
) {
	bar <- htmltools::div(
		style = list(background = fill, width = width, height = height)
	)
	chart <- htmltools::div(
		style = list(flexGrow = 1, marginLeft = "6px", background = background),
		bar
	)
	htmltools::div(
		style = list(display = "flex", alignItems = "center"), label, chart #"right"
	)
}

bar_chart_pos_neg <- function(
	label, value, max_value = 1, height = "16px",
	pos_fill = "#005ab5", neg_fill = "#dc3220"
) {
	neg_chart <- htmltools::div(style = list(flex = "1 1 0"))
	pos_chart <- htmltools::div(style = list(flex = "1 1 0"))
	width <- paste0(abs(value / max_value) * 100, "%")

	if (value < 0) {
		bar <- htmltools::div(
			style = list(
				marginLeft = "8px", background = neg_fill, width = width,
				height = height
			)
		)
		chart <- htmltools::div(
			style = list(
				display = "flex", alignItems = "center",
				justifyContent = "flex-end"
			), label, bar
		)
		neg_chart <- htmltools::tagAppendChild(neg_chart, chart)
	} else {
		bar <- htmltools::div(
			style = list(
				marginRight = "8px", background = pos_fill, width = width,
				height = height
			)
		)
		chart <- htmltools::div(
			style = list(display = "flex", alignItems = "center"), bar, label
		)
		pos_chart <- htmltools::tagAppendChild(pos_chart, chart)
	}

	htmltools::div(style = list(display = "flex"), neg_chart, pos_chart)
}


scale_check_set_default_on_fail <- function(vec, default) {
	vecf <- vec[is.finite(vec)]
	# needs at least 3 unique finite values to make a colour scale
	if(length(unique(vecf)) <= 3) {
		return(default)
	}
	# set and any infinite values to NA which can be handled by the colour scale
	if(any(!is.finite(vec))) {
		vec[is.finite(vec)] <- as.numeric(NA)
	}
	return(vec)
}

#q_scale_val <- scale_check_set_default_on_fail(data$qvalues, seq(0, 1, 0.2))

GO_GSEA_RT <- function(data) {

  if(nrow(data) == 0) {
    return(list("No results"))
  }

  q_scale_val <- data$qvalues[is.finite(data$qvalues)]
  if(length(q_scale_val) <= 1) {q_scale_val <- seq(0, 1, 0.2)}
  qvaluesf <- colourScaleR::universal_colour_scaler(
    q_scale_val,
    type = "scico", palette = "hawaii", mode = "closure", direction = 1,
    n_breaks = 9
  )

  padj_scale_val <- data$p.adjust[is.finite(data$p.adjust)]
  if(length(padj_scale_val) <= 1) {padj_scale_val <- seq(0, 1, 0.2)}
  p.adjustf <- colourScaleR::universal_colour_scaler(
    padj_scale_val,
    type = "scico", palette = "hawaii", mode = "closure", direction = 1,
    n_breaks = 9
  )

  p_scale_val <- data$pvalue[is.finite(data$pvalue)]
  if(length(p_scale_val) <= 1) {p_scale_val <- seq(0, 1, 0.2)}
  pvaluef <- colourScaleR::universal_colour_scaler(
    p_scale_val,
    type = "scico", palette = "hawaii", mode = "closure", direction = 1,
    n_breaks = 9
  )

  pc_scale_val <- data$percent[is.finite(data$percent)]
  if(length(pc_scale_val) <= 1) {pc_scale_val <- seq(0, 1, 0.2)}
  percentf <- colourScaleR::universal_colour_scaler(
    pc_scale_val,
    type = "viridis", palette = "viridis", mode = "closure", direction = 1,
    n_breaks = 9
  )

  es_up_scale_val <- data$enrichmentScore[
    data$enrichmentScore >= 0 & is.finite(data$enrichmentScore)
  ]
  if(length(es_up_scale_val) <= 1) {es_up_scale_val <- seq(0, 1, 0.2)}
  ESupf <- colourScaleR::universal_colour_scaler(
    es_up_scale_val,
    type = "brewer", palette = "RdBu", mode = "closure", direction = -1, #"RdYlBu"
    n_breaks = 8, begin = 0.6, end = 1
  )

  es_dwn_scale_val <- data$enrichmentScore[
    data$enrichmentScore < 0 & is.finite(data$enrichmentScore)
  ]
  if(length(es_dwn_scale_val) <= 1) {es_dwn_scale_val <- -seq(0, 1, 0.2)}
  ESdownf <- colourScaleR::universal_colour_scaler(
    es_dwn_scale_val,
    type = "brewer", palette = "RdBu", mode = "closure", direction = -1,
    n_breaks = 9, begin = 0, end = 0.4
  )

  nes_up_scale_val <- data$NES[data$NES >= 0 & is.finite(data$NES)]
  if(length(nes_up_scale_val) >= 0) {nes_up_scale_val <- seq(0,1,0.2)}
  NESupf <- colourScaleR::universal_colour_scaler(
    nes_up_scale_val,
    type = "brewer", palette = "RdBu", mode = "closure", direction = -1, #"RdYlBu"
    n_breaks = 8, begin = 0.6, end = 1
  )

  nes_dwn_scale_val <- data$NES[data$NES < 0 & is.finite(data$NES)]
  if(length(nes_dwn_scale_val) <= 1) {nes_dwn_scale_val <- -seq(0,1,0.2)}
  NESdownf <- colourScaleR::universal_colour_scaler(
    nes_dwn_scale_val,
    type = "brewer", palette = "RdBu", mode = "closure", direction = -1,
    n_breaks = 9, begin = 0, end = 0.4
  )

  table <-
    data %>% reactable::reactable(
      details = function(index) {
        htmltools::div(
          paste(
            "Genes: ", paste0(
              strsplit(
                unlist(data[index,"core_enrichment"]),"/"
              )[[1]], collapse = ", "
            )
          )
        )
      },
      columns = list(
        percent = reactable::colDef(
          #format = reactable::colFormat(percent = TRUE, digits = 1)
          cell = function(value) {
            valuef <- sprintf("%.1f%%", value)
            # valuef <- paste0(format(value, digits = 3),"%")
            bar_chart(valuef, width = value, fill = percentf(value))
          }#,
          #align = "right"
        ),
        enrichmentScore = reactable::colDef(
          # format = reactable::colFormat(digits = 3)
          cell = function(value) {
            bar_chart_pos_neg(
              label = format(value, digits = 3),
              value = value,
              max_value = max(abs(
                data$enrichmentScore[
                  is.finite(data$enrichmentScore)
                ]
              )),
              pos_fill = ESupf(value),
              neg_fill = ESdownf(value)
            )
          }
        ),
        NES = reactable::colDef(
          # format = reactable::colFormat(digits = 3)
          cell = function(value) {
            bar_chart_pos_neg(
              label = format(value, digits = 3),
              value = value,
              max_value = max(abs(data$NES[is.finite(data$NES)])),
              pos_fill = NESupf(value),
              neg_fill = NESdownf(value)
            )
          }
        ),
        pvalue = reactable::colDef(
          #format = reactable::colFormat(digits = 3),
          cell = function(value) {
            sprintf("%.3e", value)
          },
          style = function(value) {
            colour <- pvaluef(value)
            list(background = colour)
          }
          # bar_style(width = -log10(value))
        ),
        p.adjust = reactable::colDef(
          cell = function(value) {
            sprintf("%.3e", value)
          },
          style = function(value) {
            colour <- p.adjustf(value)
            list(background = colour)
          }
        ),
        qvalues = reactable::colDef(
          cell = function(value) {
            sprintf("%.3e", value)
          },
          style = function(value) {
            colour <- qvaluesf(value)
            list(background = colour)
          }
        ),
        core_enrichment = reactable::colDef(show = FALSE)
      ),
      showSortable = TRUE,
      searchable = TRUE, pagination = FALSE, highlight = TRUE,
      height = 960,
      filterable = TRUE
    )
  list(table)
}

GO_ORA_RT <- function(data, height = 900) {

	if(is.null(data)) {
		return(list("No results"))
	}

	if(nrow(data) == 0) {
		return(list("No results"))
	}

	q_scale_val <- scale_check_set_default_on_fail(data$qvalue, seq(0, 1, 0.2))
	qvaluef <- colourScaleR::universal_colour_scaler(
		q_scale_val,
		type = "scico", palette = "hawaii", mode = "closure", direction = 1,
		n_breaks = 9
	)

	padj_scale_val <- scale_check_set_default_on_fail(data$p.adjust, seq(0, 1, 0.2))
	p.adjustf <- colourScaleR::universal_colour_scaler(
		padj_scale_val,
		type = "scico", palette = "hawaii", mode = "closure", direction = 1,
		n_breaks = 9
	)

	p_scale_val <- scale_check_set_default_on_fail(data$pvalue, seq(0, 1, 0.2))
	pvaluef <- colourScaleR::universal_colour_scaler(
		p_scale_val,
		type = "scico", palette = "hawaii", mode = "closure", direction = 1,
		n_breaks = 9
	)

	pc_scale_val <- scale_check_set_default_on_fail(data$percent, seq(0, 1, 0.2))
	percentf <- colourScaleR::universal_colour_scaler(
		pc_scale_val,
		type = "viridis", palette = "viridis", mode = "closure", direction = 1,
		n_breaks = 9
	)

	table <-
		data %>% reactable::reactable(
			details = function(index) {
				htmltools::div(
					paste(
						"Genes: ", paste0(
							strsplit(
								unlist(data[index,"geneID"]),"/"
							)[[1]], collapse = ", "
						)
					)
				)
			},
			columns = list(
				percent = reactable::colDef(
					#format = reactable::colFormat(percent = TRUE, digits = 1)
					cell = function(value) {
						valuef <- sprintf("%.1f%%", value)
						# valuef <- paste0(format(value, digits = 3),"%")
						bar_chart(valuef, width = value, fill = percentf(value))
					}#,
					#align = "right"
				),
				pvalue = reactable::colDef(
					#format = reactable::colFormat(digits = 3),
					cell = function(value) {
						sprintf("%.3e", value)
					},
					style = function(value) {
						colour <- pvaluef(value)
						list(
							background = colour,
							color = "#FFFFFF"
						)
					}
					# bar_style(width = -log10(value))
				),
				p.adjust = reactable::colDef(
					cell = function(value) {
						sprintf("%.3e", value)
					},
					style = function(value) {
						colour <- p.adjustf(value)
						list(
							background = colour,
							color = "#FFFFFF"
						)
					}
				),
				qvalue = reactable::colDef(
					cell = function(value) {
						sprintf("%.3e", value)
					},
					style = function(value) {
						colour <- qvaluef(value)
						list(
							background = colour,
							color = "#FFFFFF"
						)
					}
				),
				GeneRatio = reactable::colDef(
					cell = function(value) {
						sprintf("%.3e", value)
					}#,
					# style = function(value) {
					# 	colour <- qvaluef(value)
					# 	list(background = colour)
					# }
				),
				BgRatio = reactable::colDef(
					cell = function(value) {
						sprintf("%.3e", value)
					}#,
					# style = function(value) {
					# 	colour <- qvaluef(value)
					# 	list(background = colour)
					# }
				),
				geneID = reactable::colDef(show = FALSE)
			),
			showSortable = TRUE,
			searchable = TRUE, pagination = FALSE, highlight = TRUE,
			height = height,
			filterable = TRUE
		)
	list(table)
}

write_FE <- function(data, name, type, out_dir) {
	fs::dir_create(out_dir)
	path <- fs::path_join(
		c(out_dir, paste0(gsub("[^\\w-]+", "_", name, perl = TRUE),".csv"))
	)
	#print(is.data.frame(GSEA_res_table_prep(data[[1]]@result)))
	#print(data[[1]]@result %>% as_tibble())
	dataf <- switch(
		type,
		"GSEA" = {
			# df <- data[[1]]@result
			# if(nrow(df) == 0) {
			# 	df <- dplyr::bind_cols(
			# 		df, tibble::tibble(
			# 			rank = integer(), leading_edge = character(),
			# 			core_enrichment = character()
			# 		)
			# 	)
			# }
			# df %>% as.data.frame() %>% GSEA_res_table_prep()
			data %>% GSEA_res_table_prep()
		},
		#"ORA" = data[[1]]@result %>% as.data.frame() %>% ORA_res_table_prep()
		"ORA" = data %>% ORA_res_table_prep()
	)
	# print(dataf)
	readr::write_csv(dataf, path)
	path
}


get_pair_names <- function(pairs, names_exp) {
	purrr::map(pairs, ~{
		paste0(
			gsub(
				"[^\\w-]+", "_",
				gsub(", ", "-", gsub("-", "_", names_exp)),
				perl = TRUE
			),
			"-", .x[2], "_vs_", .x[1]
		)
	}) %>% unlist()
}
