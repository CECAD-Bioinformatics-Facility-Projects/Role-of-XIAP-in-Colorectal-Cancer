gsInfo <- ## by getAnywhere
	function (object, geneSetID) 
	{
		geneList <- object@geneList
		if (is.numeric(geneSetID)) 
			geneSetID <- object@result[geneSetID, "ID"]
		geneSet <- object@geneSets[[geneSetID]]
		exponent <- object@params[["exponent"]]
		df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
		df$ymin <- 0
		df$ymax <- 0
		pos <- df$position == 1
		h <- diff(range(df$runningScore))/20
		df$ymin[pos] <- -h
		df$ymax[pos] <- h
		df$geneList <- geneList
		df$Description <- object@result[geneSetID, "Description"]
		df$NES         <- object@result[geneSetID, "NES"] ## UG
		df$pvalue      <- object@result[geneSetID, "pvalue"] ## UG
		df$p.adjust    <- object@result[geneSetID, "p.adjust"] ## UG
		
		return(df)
	}

gseaScores <- ## by getAnywhere 
	function (geneList, geneSet, exponent = 1, fortify = FALSE) 
{
	geneSet <- intersect(geneSet, names(geneList))
	N <- length(geneList)
	Nh <- length(geneSet)
	Phit <- Pmiss <- numeric(N)
	hits <- names(geneList) %in% geneSet
	Phit[hits] <- abs(geneList[hits])^exponent
	NR <- sum(Phit)
	Phit <- cumsum(Phit/NR)
	Pmiss[!hits] <- 1/(N - Nh)
	Pmiss <- cumsum(Pmiss)
	runningES <- Phit - Pmiss
	max.ES <- max(runningES)
	min.ES <- min(runningES)
	if (abs(max.ES) > abs(min.ES)) {
		ES <- max.ES
	}
	else {
		ES <- min.ES
	}
	df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
					 position = as.integer(hits))
	if (fortify == TRUE) {
		return(df)
	}
	df$gene = names(geneList)
	res <- list(ES = ES, runningES = df)
	return(res)
	}


gseaplot2_local <- 
function (x, geneSetID, title = "", color = "green", base_size = 11, 
		  rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
		  ES_geom = "line",
		  ## parameters by UG:
		  raw_plotlist=FALSE, linewidth = 1)  
{
	ES_geom <- match.arg(ES_geom, c("line", "dot"))
	geneList <- position <- NULL
	
	if (length(geneSetID) == 1) {
		gsdata <- gsInfo(x, geneSetID)
	}
	else {
		gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
	}

	##browser()
	p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
		theme(panel.grid.major = element_line(colour = "grey92"), 
			  panel.grid.minor = element_line(colour = "grey92"), 
			  panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
		scale_x_continuous(expand = c(0, 0))
	##browser()
	
	if (ES_geom == "line") {
		es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
							  size = linewidth)
	}
	else {
		es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
							   size = linewidth, data = subset(gsdata, position == 1))
	}
	p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
								  legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
	p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
															  axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
															  plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
															  					 unit = "cm"))
	i <- 0 
	for (term in unique(gsdata$Description)) {
		idx <- which(gsdata$ymin != 0 & gsdata$Description == 
					 	term)
		gsdata[idx, "ymin"] <- i
		gsdata[idx, "ymax"] <- i + 1
		i <- i + 1
	}
	
	## UG start -------------------------------------------------------------
	if(raw_plotlist) {
		## ---  "gsdata" is a rowwise stack of gene observations per term --- 
		## For each given panel (= barcode of one term), the x coordinate 
		## gsdata$x runs from 1 to the number of genes to be displayed 
		## (gene present in current term = bar,  absent=blank).
		## Each panel is associated with two y coordinates y_min and y_max,
		## where y_min indicates the lower y coordinate of a bar (in the current
		## panel), y_max the upper coordinate.
		## The lower or upper y coordinate can be used to count panels:
		## In the ith panel *from below!*, y_min is constant at i-1, y_max at i.
		## Likewise, gsdata$x runs from 1 to ngenes in each individual panel,
		## which means that any x in [1,ngenes] occurs in gsdata$x exactly once 
		## per panel. This, too, can aid in identifying panels.


		## Here, I count occurrences of gsdata$x == 1
		xref <- 1 
		xpos <- gsdata$x == xref ## occurrences of the reference value
		
		# use case (1): Fix the left end of the panel descriptions (to be
		#               set later)
		gsdata$x_text <- ifelse(xpos,xref,NA)

		# use case (2): Initialize a counter i with 0.
		#               Whenever the reference value xref (=1) occurs,
		#               a new panel starts. Increment i by 1  and output the 
		#               result as  the height = y-position of the current panel 
		#               description. Subtracting 0.5 shifts the text more to 
		#               the middle of each panel.
	
		i <- 0 
		gsdata$y_text <- 
			sapply(xpos, 
				   function(xp) if(xp) {
				   				  i <<- i + 1
				   				  i
				                } else {
				               	  NA
				                }
				   ) -0.5 
		
	}
	## NOTE: gsdata (both the original columns and the text coords added above)
	##       counts from y=0 upwards, i.e., the first = MOST significant term
	##       is displayed at the BOTTOM of the stacked bar charts, while the
	##       LEAST significant term is at the TOP.
	## REMEDY: In order to not mess with coordinates already stored in gsdata
	##       by the original function code, I simply reverse the entire gsdata
	##       table and invert the ymin, ymax, and y_text coordinates, such
	##       that the bar charts and their text annotation now run from the 
	##       most significant term at the top to the least significant term
	##       at the bottom.
	max_y <- max(gsdata$ymax)
	gsdata_rev <- gsdata[nrow(gsdata):1,] 
	gsdata_rev$ymin <- max_y - gsdata_rev$ymin + 1
	gsdata_rev$ymax <- max_y - gsdata_rev$ymax + 1
	gsdata_rev$y_text <- max_y - gsdata_rev$y_text + 1

	## UG end
	
	p2 <- ggplot(gsdata_rev, 
				 aes_(x = ~x)
		)+ 
		geom_linerange(aes_(ymin = ~ymin, ## original
							ymax = ~ymax, 
							color = ~Description
					   )
	    ) + 
		geom_label( ## UG
		 	aes(x = x_text,
		 		y = y_text,
		 		label = paste0("  ", ##"\t",
		 					   Description,
		 					   ## new May 8, 2026:
		 					   " [ adj. p = ",
		 					   formatC(p.adjust, 
		 					           format = "e", digits = 1), 
		 					   " ]"
		 					   ) ,
		 		hjust = "left"
		 	),
		 	size=3,
		 	alpha=0.75, #0.85,
		 	fill = "white",
		 	label.size=0.05 # border
		 		
		) +

		xlab(NULL) + ylab(NULL) + 
		theme_classic(base_size) + theme(legend.position = "none", 
										 plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
										 axis.ticks = element_blank(), 
										 axis.text = element_blank(), 
										 axis.line.x = element_blank()
										 ) + 
		scale_x_continuous(expand = c(0, 0)) + 
		scale_y_continuous(expand = c(0, 0))


	if (length(geneSetID) == 1) {
		v <- seq(1, sum(gsdata$position), length.out = 9)
		inv <- findInterval(rev(cumsum(gsdata$position)), v)
		if (min(inv) == 0) 
			inv <- inv + 1
		col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
		ymin <- min(p2$data$ymin)
		yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
		xmin <- which(!duplicated(inv))
		xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
		d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
						xmax = xmax, col = col[unique(inv)])
		p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
								  ymin = ~ymin, ymax = ~ymax, fill = ~I(col)
						     ), 
							 data = d, 
							 alpha = 0.9, inherit.aes = FALSE)
	}
	df2 <- p$data
	df2$y <- p$data$geneList[df2$x]
	p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
											   y = ~y, yend = 0), 
							  color = "grey",
							  size=2
							  ) 
	                          

	p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + ## original

		theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
								   l = 0.2, unit = "cm"))
	if (!is.null(title) && !is.na(title) && title != "") 
		p.res <- p.res + ggtitle(title)
	if (length(color) == length(geneSetID)) {
		p.res <- p.res + scale_color_manual(values = color)
		if (length(color) == 1) {
			p.res <- p.res + theme(legend.position = "none")
			p2 <- p2 + scale_color_manual(values = "black")
		}
		else {
			p2 <- p2 + scale_color_manual(values = color)
		}
	}
	if (pvalue_table) {
		pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
		rownames(pd) <- pd$Description
		pd <- pd[, -1]
		pd <- round(pd, 4)
		tp <- tableGrob2(pd, p.res)
		p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
																			 xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
																			 													0.95), ymin = quantile(p.res$data$runningScore, 
																			 																		   0.75), ymax = quantile(p.res$data$runningScore, 
																			 																		   					   0.9))
	}
	plotlist <- list(p.res, p2, p.pos)[subplots]
	
	if(raw_plotlist) { ## UG
		return(plotlist)
	}
	
	n <- length(plotlist)
	plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
										   axis.ticks.x = element_line(), axis.text.x = element_text())
	if (length(subplots) == 1) 
		return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
														  r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
	if (length(rel_heights) > length(subplots)) 
		rel_heights <- rel_heights[subplots]
	aplot::plot_list(gglist = plotlist, ncol = 1, heights = rel_heights)
}


