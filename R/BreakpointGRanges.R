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
#'@export
partner <- function(gr) {
	assertthat::assert_that(all(gr$partner %in% names(gr)))
	return(gr[gr$partner,])
}

#' Extracts the breakpoint sequence
#' 
#' @param gr breakpoint GRanges
#' @param ref Reference BSgenome
#' @param anchoredBases Number of bases leading into breakpoint to extract
#' @param remoteBases Number of bases from other side of breakpoint to extract
#' @param position breakpoint position to use. In the presence of homology
#'     or imprecise calls, multiple breakpoint positions are possible.
#'     Valid values are "start", "end", "middle".
#' @export
breakpointSequence <- function(gr, ref, anchoredBases=50, remoteBases=anchoredBases, position="middle") {
}
#' Returns the reference sequence around the breakpoint position
#' 
#' @param gr breakpoint GRanges
#' @param ref Reference BSgenome
#' @param anchoredBases Number of bases leading into breakpoint to extract
#' @param followingBases Number of reference bases past breakpoint to extract
#' @param position breakpoint position to use. In the presence of homology
#'     or imprecise calls, multiple breakpoint positions are possible.
#'     Valid values are "start", "end", "middle".
#' @export
referenceSequence <- function(gr, ref, anchoredBases=50, followingBases=anchoredBases, position="middle") {
}

#'@export
referenceHomology <- function(gr, ref,
		flankSize=100,
		margin=5,
		match = 2, mismatch = -6, gapOpening = 5, gapExtension = 3, # bwa
		#match = 1, mismatch = -4, gapOpening = 6, gapExtension = 1, # bowtie2
		...) {
	varseq <- breakpointSequence(gr, ref, anchoredBases=flankSize)
	refseq <- referenceSequence(gr, ref, anchoredBases=flankSize, followingBases=nchar(varseq) - 2 * flankSize + margin)
	aln <- pairwiseAlignment(varseq, refseq, type="global-local",
 		substitutionMatrix=nucleotideSubstitutionMatrix(match, mismatch, FALSE, "DNA"),
 		gapOpening=gapOpening, gapExtension=gapExtension)
	partnerIndex <- match(gr$partner, names(gr))
	homlen <- nchar(aln) - deletion(nindel(aln))[,2]- insertion(nindel(aln))[,2]
	bphomlen <- homlen + homlen[partnerIndex]
	bpscore <- score(aln) + score(aln)[partnerIndex] - 2 * flankSize * match
	return(list(homlen=bphomlen, score=bpscore))
}

blastHomology <- function(gr, ...) {
	requireNamespace("rBLAST", quietly=FALSE)
}