#' Converts a breakpoint GRanges to a BEDPE-formatted data frame
#' @param gr GRanges object
#'
#' @return BEDPE-formatted data frame
#' @export
breakpointgr2bedpe <- function(gr){
	bedpe <- data.frame(
		chrom1=seqnames(gr),
		start1=start(gr) - 1,
		end1=end(gr),
		chrom2=seqnames(partner(gr)),
		start2=start(partner(gr)) - 1,
		end2=end(partner(gr)),
		name=names(gr),
		partner.name=names(partner(gr)),
		score=gr$QUAL,
		strand1=strand(gr),
		strand2=strand(partner(gr))
	)
	l <- mapply(c, bedpe$name, bedpe$partner.name, SIMPLIFY = FALSE)
	ll <- mapply(sort, l, SIMPLIFY = FALSE)
	bedpe <- bedpe %>% dplyr::mutate(pair = as.character(ll)) %>%
		dplyr::distinct(pair,.keep_all = TRUE) %>% dplyr::select(-c('partner.name','pair'))
	#bedpe <- bedpe[stringr::str_detect(bedpe$name, "gridss.+o"),]
	return(bedpe)
}

#' @param bedpe BEDPE-formatted dataframe
#' @param ... plotting parameters for circlize::circos.initializeWithIdeogram() and circlize::circos.genomicLink()
#' @return circos plot with connected breakpoints
bedpe2circos <- function(bedpe,
						 cytoband = system.file(
						 	package = "circlize","extdata", "cytoBand.txt"),
						 species = NULL, sort.chr = TRUE,
						 chromosome.index = NULL, major.by = NULL,
						 plotType = c("ideogram", "axis", "labels"),
						 track.height = NULL,
						 ideogram.height = circlize::convert_height(2, "mm"),
						 rou = circlize:::get_most_inside_radius(),
						 rou1 = rou, rou2 = rou,
						 col = "black", lwd = par("lwd"),
						 lty = par("lty"), border = col, ...){
	`%<a-%` <- pryr::`%<a-%`
	bed1<-bedpe[,1:3]
	bed2<-bedpe[,4:6]
	circos.pryr %<a-% {
		circlize::circos.initializeWithIdeogram(
			cytoband = cytoband,
				species = species, sort.chr = sort.chr,
				chromosome.index = chromosome.index, major.by = major.by,
				plotType = plotType,
				track.height = track.height, ideogram.height = ideogram.height, ...)
		circlize::circos.genomicLink(bed1, bed2, rou = rou,
									 rou1 = rou1, rou2 = rou2,
									 col = col, lwd = lwd,
									 lty = lty, border = border, ...)

	}
	circlize::circos.clear()
	return(circos.pryr)
}

