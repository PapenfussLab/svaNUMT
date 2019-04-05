#' GRanges representing the breakend coordinates of
#' structural variants
#' #@export
#setClass("BreakpointGRanges", contains="GRanges")

#' Partner breakend for each breakend.
#'
#' @details
#' All breakends must have their partner breakend included
#' in the GRanges.
#'
#' @param gr GRanges object of SV breakends
#' @return A GRanges object in which each entry is the partner breakend of
#' those in the input object.
#' @examples
#' #reading in a VCF file as \code{vcf}
#' vcf.file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19") vcf.file <- system.file("extdata", "g
#' #parsing \code{vcf} to GRanges object \code{gr}
#' gr <- breakpointRanges(vcf)
#' #output partner breakend of each breakend in \code{gr}
#' partner(gr)
#'@export
partner <- function(gr) {
	assertthat::assert_that(all(gr$partner %in% names(gr)))
	return(gr[gr$partner,])
}

#' Finding overlapping breakpoints between two breakpoint sets
#'
#' @details
#' \code{findBreakpointOverlaps()} is an efficient adaptation of \code{findOverlaps-methods()}
#' for breakend ranges. It searches for overlaps between breakpoint objects, and return a
#' matrix including index of overlapping ranges as well as error stats.
#' All breakends must have their partner breakend included in the \code{partner}
#' field. A valid overlap requires that breakends on boths sides meets the overlapping
#' requirements.
#'
#' See GenomicRanges::findOverlaps-methods for details of overlap calculation.
#'
#' @param query,subject Both of the input objects should be GRanges objects.
#' Unlike \code{findOverlaps()}, \code{subject} cannot be ommitted. Each breakpoint
#' must be accompanied with a partner breakend, which is also in the GRanges, with the
#' partner's id recorded in the \code{partner} field.
#' See GenomicRanges::findOverlaps-methods for details.
#' @param maxgap,minoverlap Valid overlapping thresholds of a maximum gap and a minimum
#' overlapping positions between breakend intervals. Both should be scalar integers. maxgap
#' allows non-negative values, and minoverlap allows positive values.
#' See GenomicRanges::findOverlaps-methods for details.
#' @param ignore.strand Default value is FALSE. strand inforamtion is ignored when set to
#' TRUE.
#' See GenomicRanges::findOverlaps-methods for details.
#' @param sizemargin Error margin in allowable size to prevent matching of events
#' of different sizes, e.g. a 200bp event matching a 1bp event when maxgap is
#' set to 200.
#' @param restrictMarginToSizeMultiple Size restriction multiplier on event size.
#' The default value of 0.5 requires that the breakpoint positions can be off by
#' at maximum, half the event size. This ensures that small deletion do actually
#' overlap at least one base pair.
#' @examples
#' #reading in VCF files
#' query.file <- system.file("extdata", "gridss-na12878.vcf", package = "StructuralVariantAnnotation")
#' subject.file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
#' query.vcf <- VariantAnnotation::readVcf(query.file, "hg19")
#' subject.vcf <- VariantAnnotation::readVcf(subject.file, "hg19")
#' #parsing vcfs to GRanges objects
#' query.gr <- breakpointRanges(query.vcf)
#' subject.gr <- breakpointRanges(subject.vcf)
#' #find overlapping breakpoint intervals
#' findBreakpointOverlaps(query.gr, subject.gr)
#' findBreakpointOverlaps(query.gr, subject.gr, ignore.strand=TRUE)
#' findBreakpointOverlaps(query.gr, subject.gr, maxgap=100, minoverlap=10, sizemargin=0.5)
#' @return A dataframe containing index and error stats of overlapping breakpoints.
#'@export
findBreakpointOverlaps <- function(query, subject, maxgap=-1L, minoverlap=0L, ignore.strand=FALSE, sizemargin=0.25, restrictMarginToSizeMultiple=0.5) {
	hits <- dplyr::bind_rows(
		as.data.frame(findOverlaps(query, subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand), row.names=NULL),
		as.data.frame(findOverlaps(partner(query), partner(subject), maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand), row.names=NULL))
	#as.data.frame(findOverlaps(partner(query), partner(subject), maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand), row.names=NULL))

	# we now want to do:
	# hits <- hits[duplicated(hits),] # both breakends match
	# but for large hit sets (such as focal false positive loci) we run out of memory (>32GB)
	# instead, we sort then check that we match the next record
	hits = hits %>% dplyr::arrange(queryHits, subjectHits) %>%
		dplyr::filter(!is.na(dplyr::lead(.$queryHits)) & !is.na(dplyr::lead(.$subjectHits)) & dplyr::lead(.$queryHits) == .$queryHits & dplyr::lead(.$subjectHits) == .$subjectHits)
	if (!is.null(sizemargin) && !is.na(sizemargin)) {
		# take into account confidence intervals when calculating event size
		callwidth <- .distance(query, partner(query))
		truthwidth <- .distance(subject, partner(subject))
		callsize <- callwidth + (query$insLen %na% 0)
		truthsize <- truthwidth + (subject$insLen %na% 0)
		hits$sizeerror <- .distance(
			IRanges(start=callsize$min[hits$queryHits], end=callsize$max[hits$queryHits]),
			IRanges(start=truthsize$min[hits$subjectHits], end=truthsize$max[hits$subjectHits])
			)$min
		# event sizes must be within sizemargin
		hits <- hits[hits$sizeerror - 1 < sizemargin * pmin(callsize$max[hits$queryHits], truthsize$max[hits$subjectHits]),]
		# further restrict breakpoint positions for small events
		hits$localbperror <- .distance(query[hits$queryHits], subject[hits$subjectHits])$min
		hits$remotebperror <- .distance(partner(query)[hits$queryHits], partner(subject)[hits$subjectHits])$min
		if (!is.null(restrictMarginToSizeMultiple)) {
			allowablePositionError <- (pmin(callsize$max[hits$queryHits], truthsize$max[hits$subjectHits]) * restrictMarginToSizeMultiple + 1)
			hits <- hits[hits$localbperror <= allowablePositionError & hits$remotebperror <= allowablePositionError, ]
		}
	}
	return(hits)
}
.distance <- function(r1, r2) {
	return(data.frame(
		min=pmax(0, pmax(start(r1), start(r2)) - pmin(end(r1), end(r2))),
		max=pmax(end(r2) - start(r1), end(r1) - start(r2))))
}
#' Counting overlapping breakpoints between two breakpoint sets
#'
#' @details
#' \code{countBreakpointOverlaps()} returns the number of overlaps between breakpoint
#' objects, based on the output of \code{findBreakpointOverlaps()}.
#' See GenomicRanges::countOverlaps-methods
#' @param querygr,subjectgr,maxgap,minoverlap,ignore.strand,sizemargin,restrictMarginToSizeMultiple
#' See \code{findBreakpointOverlaps()}.
#' @param countOnlyBest Default value set to FALSE. When set to TRUE, the result count
#' each subject breakpoint as overlaping only the best overlapping query breakpoint.
#' @param breakpointScoreColumn Query column defining a score for determining which query breakpoint
#' is considered the best when countOnlyBest=TRUE.
#' @examples
#' #reading in VCF files
#' query.file <- system.file("extdata", "gridss-na12878.vcf", package = "StructuralVariantAnnotation")
#' subject.file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
#' query.vcf <- VariantAnnotation::readVcf(query.file, "hg19")
#' subject.vcf <- VariantAnnotation::readVcf(subject.file, "hg19")
#' #parsing vcfs to GRanges objects
#' query.gr <- breakpointRanges(query.vcf)
#' subject.gr <- breakpointRanges(subject.vcf)
#' #count overlapping breakpoint intervals
#' countBreakpointOverlaps(query.gr, subject.gr)
#' countBreakpointOverlaps(query.gr, subject.gr, maxgap=100)
#' countBreakpointOverlaps(query.gr, subject.gr, maxgap=100, ignore.strand=TRUE, countOnlyBest=TRUE)
#' @return An integer vector containing the tabulated query overlap hits.
#' @export
countBreakpointOverlaps <- function(querygr, subjectgr, countOnlyBest=FALSE,
									breakpointScoreColumn = "QUAL", maxgap=-1L,
									minoverlap=0L, ignore.strand=FALSE, sizemargin=0.25,
									restrictMarginToSizeMultiple=0.5) {
	hitscounts <- rep(0, length(querygr))
	hits <- findBreakpointOverlaps(querygr, subjectgr, maxgap, minoverlap, ignore.strand, sizemargin=sizemargin, restrictMarginToSizeMultiple=restrictMarginToSizeMultiple)
	if (!countOnlyBest) {
		hits <- hits %>%
	      dplyr::group_by(queryHits) %>%
	      dplyr::summarise(n=dplyr::n())
	} else {
		# assign supporting evidence to the call with the highest QUAL
		hits$QUAL <- mcols(querygr)[[breakpointScoreColumn]][hits$queryHits]
	    hits <- hits %>%
	      dplyr::arrange(desc(QUAL), queryHits) %>%
	      dplyr::distinct(subjectHits, .keep_all=TRUE) %>%
	      dplyr::group_by(queryHits) %>%
	      dplyr::summarise(n=n())
	}
    hitscounts[hits$queryHits] <- hits$n
    return(hitscounts)
}

#' Loading a breakpoint GRanges from a BEDPE file
#' @details
#' A BEDPE file is taken as an input and converted to a breakpoint GRanges object. The file is first loaded using
#' rtracklayer::import as a Pairs object, then converted to GRanges.
#' @param file A BEDPE file. See \url{https://bedtools.readthedocs.io/en/latest/content/general-usage.html} for details.
#' @param placeholderName Prefix to use to ensure each entry has a unique ID.
#' @examples
#' bedpe.file <- system.file("extdata", "gridss.bedpe", package = "StructuralVariantAnnotation")
#' bedpe2breakpointgr(bedpe.file)
#' bedpe2breakpointgr(bedpe.file, placeholderName='gridss')
#' @return breakpoint GRanges object
#' @export
bedpe2breakpointgr <- function(file, placeholderName="bedpe") {
	return(pairs2breakpointgr(import(file), placeholderName))
}
#' Converts a BEDPE Pairs containing pairs of GRanges loaded using to a breakpoint GRanges object.
#' @param pairs a Pairs object consisting of two parallel genomic loci.
#' @param placeholderName prefix to use to ensure each entry has a unique ID.
#' @examples
#' bedpe.file <- system.file("extdata", "gridss.bedpe", package = "StructuralVariantAnnotation")
#' bedpe.pairs <- rtracklayer::import(bedpe.file)
#' pairs2breakpointgr(bedpe.pairs)
#' pairs2breakpointgr(bedpe.pairs, placeholderName='gridss')
#' @return breakpoint GRanges object
#' @export
pairs2breakpointgr <- function(pairs, placeholderName="bedpe") {
	n <- names(pairs)
	if (is.null(n)) {
		# BEDPE uses the "name" field
		if ("name" %in% names(mcols(pairs))) {
			n <- mcols(pairs)$name
		} else {
			n <- rep(NA_character_, length(pairs))
		}
	}
	# ensure row names are unique
	n <- ifelse(is.na(n) | n == "" | n =="." | duplicated(n), paste0(placeholderName, seq_along(n)), n)
	#
	gr <- c(S4Vectors::first(pairs), S4Vectors::second(pairs))
	names(gr) <- c(paste0(n, "_1"), paste0(n, "_2"))
	gr$partner <- c(paste0(n, "_2"), paste0(n, "_1"))
	for (col in names(mcols(pairs))) {
		if (col %in% c("name")) {
			# drop columns we have processed
		} else {
			mcols(gr)[[col]] <- mcols(pairs)[[col]]
		}
	}
	return(gr)
}
