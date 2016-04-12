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

#' Finds overlapping breakpoints by requiring that breakends on
#' boths sides overlap
#'
#' @details
#' See GenomicRanges::findOverlaps-methods for details of overlap calculation
#'
#'@export
findBreakpointOverlaps <- function(query, subject, maxgap=0L, minoverlap=1L, ignore.strand=FALSE) {
  dfhits <- rbind(as.data.frame(findOverlaps(query, subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand), row.names=NULL),
                  as.data.frame(findOverlaps(partner(query), partner(subject), maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand), row.names=NULL))
  dfhits <- dfhits[duplicated(dfhits),] # both breakends match
  row.names(dfhits) <- NULL
  return(dfhits)
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
	gr <- .constrict(gr)
	seqgr <- GRanges(seqnames=seqnames(gr), ranges=IRanges(
		start=start(gr) - ifelse(strand(gr) == "-", followingBases, anchoredBases - 1),
		end=end(gr) + ifelse(strand(gr) == "-", anchoredBases - 1, followingBases)))
	startPad <- pmax(0, 1 - start(seqgr))
	endPad <- pmax(0, end(seqgr) - seqlengths(ref)[as.character(seqnames(seqgr))])
	ranges(seqgr) <- IRanges(start=start(seqgr) + startPad, end=end(seqgr) - endPad)
	seq <- Biostrings::getSeq(ref, seqgr)
	seq <- paste0(stringr::str_pad("", startPad, pad="N"), as.character(seq), stringr::str_pad("", endPad, pad="N"))
	# DNAStringSet doesn't like out of bounds subsetting
	seq <- ifelse(strand(gr) == "-", as.character(Biostrings::reverseComplement(DNAStringSet(seq))), seq)
	return(seq)
}
#' constrict
.constrict <- function(gr, ref=NULL,position="middle") {
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

	} else {
		stop(paste("Unrecognised position", position))
	}
	if (!is.null(ref)) {
		ranges(gr) <- IRanges(start=pmin(pmax(1, start(gr)), seqlengths(ref)[as.character(seqnames(gr))]), width=1)
	}
	return(gr)
}

#' Calculates the length of inexact homology between the breakpoint sequence
#' and the reference
#'
#' @param gr breakpoint GRanges
#' @param ref Reference BSgenome
#' @param anchorLength Number of bases to consider for homology
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
		anchorLength=300,
		margin=5,
		match=2, mismatch=-6, gapOpening=5, gapExtension=3 # bwa
		#match = 1, mismatch = -4, gapOpening = 6, gapExtension = 1, # bowtie2
		) {
	# shrink anchor for small events to prevent spanning alignment
	aLength <- pmin(anchorLength, abs(gr$svLen) + 1) %na% anchorLength
	anchorSeq <- referenceSequence(gr, ref, aLength, 0)
	anchorSeq <- sub(".*N", "", anchorSeq)
	# shrink anchor with Ns
	aLength <- nchar(anchorSeq)
	varseq <- breakpointSequence(gr, ref, aLength)
	varseq <- sub("N.*", "", varseq)
	bpLength <- nchar(varseq) - aLength
	nonbpseq <- referenceSequence(gr, ref, 0, bpLength + margin)
	nonbpseq <- sub("N.*", "", nonbpseq)
	refseq <- paste0(anchorSeq, nonbpseq)

	partnerIndex <- match(gr$partner, names(gr))

	if (all(refseq=="") && all(varseq=="")) {
		# Workaround of Biostrings::pairwiseAlignment bug
		return(data.frame(
			exacthomlen=rep(NA, length(gr)),
			inexacthomlen=rep(NA, length(gr)),
			inexactscore=rep(NA, length(gr))))
	}

	aln <- Biostrings::pairwiseAlignment(varseq, refseq, type="local",
 		substitutionMatrix=nucleotideSubstitutionMatrix(match, mismatch, FALSE, "DNA"),
 		gapOpening=gapOpening, gapExtension=gapExtension, scoreOnly=FALSE)
	ihomlen <- Biostrings::nchar(aln) - aLength - deletion(nindel(aln))[,2] - insertion(nindel(aln))[,2]
	ibphomlen <- ihomlen + ihomlen[partnerIndex]
	ibpscore <- score(aln) + score(aln)[partnerIndex] - 2 * aLength * match

	# TODO: replace this with an efficient longest common substring function
	# instead of S/W with a massive mismatch/gap penalty
	penalty <- anchorLength * match
	matchLength <- Biostrings::pairwiseAlignment(varseq, refseq, type="local",
 		substitutionMatrix=nucleotideSubstitutionMatrix(match, -penalty, FALSE, "DNA"),
 		gapOpening=penalty, gapExtension=0, scoreOnly=TRUE) / match
	ehomlen <- matchLength - aLength
	ebphomlen <- ehomlen + ehomlen[partnerIndex]

	ebphomlen[aLength == 0] <- NA
	ibphomlen[aLength == 0] <- NA
	ibpscore[aLength == 0] <- NA
	return(data.frame(
		exacthomlen=ebphomlen,
		inexacthomlen=ibphomlen,
		inexactscore=ibpscore))
}

#' Identifies breakpoint sequences with signficant homology to BLAST database
#' sequences. Apparent breakpoints containing such sequence are better explained
#' by the sequence from the BLAST database such as by alternate assemblies.
#'
#' @details
#' See https://github.com/mhahsler/rBLAST for rBLAST package installation
#' instructions
#' Download and install the package from AppVeyor or install via install_github("mhahsler/rBLAST") (requires the R package devtools)
blastHomology <- function(gr, ref, db, anchorLength=150) {
	requireNamespace("rBLAST", quietly=FALSE)
	blastseq <- DNAStringSet(breakpointSequence(gr, ref, anchorLength))
	bl <- rBLAST::blast(db=db)
	cl <- predict(bl, blastseq)
	cl$index <- as.integer(substring(cl$QueryID, 7))
	cl$leftOverlap <- anchorLength - cl$Q.start + 1
	cl$rightOverlap <- cl$Q.end - (nchar(blastseq) - anchorLength)
	cl$minOverlap <- pmin(cl$leftOverlap, cl$rightOverlap)
	return(cl)
}

