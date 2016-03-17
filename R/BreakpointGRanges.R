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

#' Extracts the breakpoint sequence.
#' 
#' @details
#' The sequence is the sequenced traversed from the reference anchor bases
#' to the breakpoint. For backward (-) breakpoints, this corresponds to the
#' reverse compliment of the reference sequence bases.
#' 
#' @param gr breakpoint GRanges
#' @param ref Reference BSgenome
#' @param anchoredBases Number of bases leading into breakpoint to extract
#' @param remoteBases Number of bases from other side of breakpoint to extract
#' @export
breakpointSequence <- function(gr, ref, anchoredBases, remoteBases=anchoredBases) {
	localSeq <- referenceSequence(gr, ref, anchoredBases, 0)
	insSeq <- ifelse(strand(gr) == "-",
		as.character(Biostrings::reverseComplement(DNAStringSet(gr$insSeq %na% ""))),
		gr$insSeq %na% "")
	remoteSeq <- as.character(Biostrings::reverseComplement(DNAStringSet(
		referenceSequence(partner(gr), ref, remoteBases, 0))))
	return(paste0(localSeq, insSeq, remoteSeq))
}
#' Returns the reference sequence around the breakpoint position
#' 
#' @details
#' The sequence is the sequenced traversed from the reference anchor bases
#' to the breakpoint. For backward (-) breakpoints, this corresponds to the
#' reverse compliment of the reference sequence bases.
#' 
#' @param gr breakpoint GRanges
#' @param ref Reference BSgenome
#' @param anchoredBases Number of bases leading into breakpoint to extract
#' @param followingBases Number of reference bases past breakpoint to extract
#' @export
referenceSequence <- function(gr, ref, anchoredBases, followingBases=anchoredBases) {
	assertthat::assert_that(is(gr, "GRanges"))
	assertthat::assert_that(is(ref, "BSgenome"))
	seqgr <- GRanges(seqnames=seqnames(gr), ranges=IRanges(
		start=start(gr) - ifelse(strand(gr) == "-", followingBases, anchoredBases - 1),
		end=end(gr) + ifelse(strand(gr) == "-", anchoredBases - 1, followingBases)))
	seq <- Biostrings::getSeq(ref, seqgr)
	# DNAStringSet doesn't like out of bounds subsetting
	seq <- ifelse(strand(gr) == "-", as.character(Biostrings::reverseComplement(seq)), as.character(seq))
	return(seq)
}
#' constrict
.constrict <- function(gr, position="middle") {
	isLower <- start(gr) < start(partner(gr))
	# Want to call a valid breakpoint
	#  123 456
	#
	#  =>   <= + -
	#  >   <== f f
	#
	#  =>  =>  + +
	#  >   ==> f c
	roundDown <- isLower | strand(gr) == "-"
	if (position == "middle") {
		pos <- (start(gr) + end(gr)) / 2
		ranges(gr) <- IRanges(
			start=ifelse(roundDown,floor(pos), ceiling(pos)),
			width=1, names=names(gr))
		return(gr)
	}
	stop(paste("Unrecognised position", position))
}


#' Calculates the length of inexact homology between the breakpoint sequence
#' and the reference
#' 
#' @param gr breakpoint GRanges
#' @param ref Reference BSgenome
#' @param flankSize Number of bases to consider for homology
#' @param margin Number of additional reference bases include. This allows
#'		for inexact homology to be detected even in the presence of indels.
#' @param match alignment 
#' @param mismatch see Biostrings::pairwiseAlignment
#' @param gapOpening see Biostrings::pairwiseAlignment
#' @param gapExtension see Biostrings::pairwiseAlignment
#' @param match see Biostrings::pairwiseAlignment
#'     
#'@export
referenceHomology <- function(gr, ref,
		flankSize=200,
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

