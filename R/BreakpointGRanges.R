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
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
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
#' @param ignore.strand Default value is FALSE. strand information is ignored when set to
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
#' findBreakpointOverlaps(query.gr, subject.gr, maxgap=100, sizemargin=0.5)
#' @return A dataframe containing index and error stats of overlapping breakpoints.
#'@export
findBreakpointOverlaps <- function(query, subject, maxgap=-1L, minoverlap=0L, ignore.strand=FALSE, sizemargin=NULL, restrictMarginToSizeMultiple=NULL) {
	localhits = findOverlaps(query, subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand)
	remotehits = findOverlaps(partner(query), partner(subject), maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand)
	## duplicated() version:
	#hits = Hits(c(queryHits(localhits), queryHits(remotehits)), c(subjectHits(localhits), subjectHits(remotehits)), nLnode=nLnode(localhits), nRnode=nRnode(localhits), sort.by.query=TRUE)
	#hits = hits[duplicated(hits)]
	
	## intersect() version:
	hits = BiocGenerics::intersect(localhits, remotehits)
	
	## dplyr() version:
	#hits <- dplyr::bind_rows(
	#	as.data.frame(localhits, row.names=NULL),
	#	as.data.frame(remotehits, row.names=NULL))
	#hits = hits %>% dplyr::arrange(queryHits, subjectHits) %>%
	#	dplyr::filter(!is.na(dplyr::lead(.$queryHits)) & !is.na(dplyr::lead(.$subjectHits)) & dplyr::lead(.$queryHits) == .$queryHits & dplyr::lead(.$subjectHits) == .$subjectHits)
	
	## dplyr() exploiting the sorted nature of the findOverlaps():
	#hits = Hits(c(queryHits(localhits), queryHits(remotehits)), c(subjectHits(localhits), subjectHits(remotehits)), nLnode=nLnode(localhits), nRnode=nRnode(localhits), sort.by.query=TRUE)
	#queryLead  = dplyr::lead(queryHits(hits))
	#querySubject  = dplyr::lead(queryHits(hits))
	#hits = hits[
	#	!is.na(queryLead) &d
	#	!is.na(querySubject) &
	#	queryLead == queryHits(hits) &
	#	querySubject == subjectHits(hits)]
	
	if (!is.null(sizemargin) && !is.na(sizemargin)) {
		# take into account confidence intervals when calculating event size
		callwidth <- .distance(query, partner(query))
		truthwidth <- .distance(subject, partner(subject))
		callsize <- callwidth + (query$insLen %na% 0)
		truthsize <- truthwidth + (subject$insLen %na% 0)
		sizeerror <- .distance(
			IRanges::IRanges(start=callsize$min[S4Vectors::queryHits(hits)], end=callsize$max[S4Vectors::queryHits(hits)]),
			IRanges::IRanges(start=truthsize$min[S4Vectors::subjectHits(hits)], end=truthsize$max[S4Vectors::subjectHits(hits)])
			)$min
		# event sizes must be within sizemargin
		hits <- hits[sizeerror - 1 < sizemargin * pmin(callsize$max[S4Vectors::queryHits(hits)], truthsize$max[S4Vectors::subjectHits(hits)]),]
		# further restrict breakpoint positions for small events
		localbperror <- .distance(query[S4Vectors::queryHits(hits)], subject[S4Vectors::subjectHits(hits)])$min
		remotebperror <- .distance(partner(query)[S4Vectors::queryHits(hits)], partner(subject)[S4Vectors::subjectHits(hits)])$min
		if (!is.null(restrictMarginToSizeMultiple)) {
			allowablePositionError <- (pmin(callsize$max[S4Vectors::queryHits(hits)], truthsize$max[S4Vectors::subjectHits(hits)]) * restrictMarginToSizeMultiple + 1)
			hits <- hits[localbperror <= allowablePositionError & remotebperror <= allowablePositionError, ]
		}
	}
	return(hits)
}
# TODO: new function to annotate a Hits object with sizeerror, localbperror, and remotebperror
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
									minoverlap=0L, ignore.strand=FALSE, sizemargin=NULL,
									restrictMarginToSizeMultiple=NULL) {
	hitscounts <- rep(0, length(querygr))
	hits <- as.data.frame(findBreakpointOverlaps(querygr, subjectgr, maxgap, minoverlap, ignore.strand, sizemargin=sizemargin, restrictMarginToSizeMultiple=restrictMarginToSizeMultiple))
	if (!countOnlyBest) {
		hits <- hits %>%
	      dplyr::group_by(.data$queryHits) %>%
	      dplyr::summarise(n=dplyr::n())
	} else {
		# assign supporting evidence to the call with the highest QUAL
		hits$QUAL <- S4Vectors::mcols(querygr)[[breakpointScoreColumn]][hits$queryHits]
	    hits <- hits %>%
	      dplyr::arrange(desc(.data$QUAL), .data$queryHits) %>%
	      dplyr::distinct(.data$subjectHits, .keep_all=TRUE) %>%
	      dplyr::group_by(.data$queryHits) %>%
	      dplyr::summarise(n=dplyr::n())
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
		if ("name" %in% names(S4Vectors::mcols(pairs))) {
			n <- S4Vectors::mcols(pairs)$name
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
	for (col in names(S4Vectors::mcols(pairs))) {
		if (col %in% c("name")) {
			# drop columns we have processed
		} else {
			S4Vectors::mcols(gr)[[col]] <- S4Vectors::mcols(pairs)[[col]]
		}
	}
	return(gr)
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
extractBreakpointSequence <- function(gr, ref, anchoredBases, remoteBases=anchoredBases) {
	localSeq <- extractReferenceSequence(gr, ref, anchoredBases, 0)
	insSeq <- ifelse(strand(gr) == "-",
					 as.character(Biostrings::reverseComplement(DNAStringSet(gr$insSeq %na% ""))),
					 gr$insSeq %na% "")
	remoteSeq <- as.character(Biostrings::reverseComplement(DNAStringSet(
		extractReferenceSequence(partner(gr), ref, remoteBases, 0))))
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
extractReferenceSequence <- function(gr, ref, anchoredBases, followingBases=anchoredBases) {
	assertthat::assert_that(is(gr, "GRanges"))
	assertthat::assert_that(is(ref, "BSgenome"))
	gr <- .constrict(gr)
	seqgr <- GRanges(seqnames=GenomeInfoDb::seqnames(gr), ranges=IRanges::IRanges(
		start=start(gr) - ifelse(strand(gr) == "-", followingBases, anchoredBases - 1),
		end=end(gr) + ifelse(strand(gr) == "-", anchoredBases - 1, followingBases)))
	startPad <- pmax(0, 1 - start(seqgr))
	endPad <- pmax(0, end(seqgr) - GenomeInfoDb::seqlengths(ref)[as.character(GenomeInfoDb::seqnames(seqgr))])
	GenomicRanges::ranges(seqgr) <- IRanges::IRanges(start=start(seqgr) + startPad, end=end(seqgr) - endPad)
	seq <- Biostrings::getSeq(ref, seqgr)
	seq <- paste0(stringr::str_pad("", startPad, pad="N"), as.character(seq), stringr::str_pad("", endPad, pad="N"))
	# DNAStringSet doesn't like out of bounds subsetting
	seq <- ifelse(strand(gr) == "-", as.character(Biostrings::reverseComplement(DNAStringSet(seq))), seq)
	return(seq)
}
#' constrict
#' @param gr GRanges object
#' @param ref reference 
#' @param position only 'middle' position is accepted.
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
		GenomicRanges::ranges(gr) <- IRanges::IRanges(
			start=ifelse(roundDown,floor(pos), ceiling(pos)),
			width=1, names=names(gr))

	} else {
		stop(paste("Unrecognised position", position))
	}
	if (!is.null(ref)) {
		GenomicRanges::ranges(gr) <- IRanges::IRanges(start=pmin(pmax(1, start(gr)), GenomeInfoDb::seqlengths(ref)[as.character(GenomeInfoDb::seqnames(gr))]), width=1)
	}
	return(gr)
}

#' Calculates the length of inexact homology between the breakpoint sequence
#' and the reference
#'
#' @param gr reakpoint GRanges
#' @param ref reference BSgenome
#' @param anchorLength Number of bases to consider for homology
#' @param margin Number of additional reference bases include. This allows
#'		for inexact homology to be detected even in the presence of indels.
#' @param mismatch see Biostrings::pairwiseAlignment
#' @param gapOpening see Biostrings::pairwiseAlignment
#' @param gapExtension see Biostrings::pairwiseAlignment
#' @param match see Biostrings::pairwiseAlignment
#'
calculateReferenceHomology <- function(gr, ref,
									   anchorLength=300,
									   margin=5,
									   match=2, mismatch=-6, gapOpening=5, gapExtension=3 # bwa
									   #match = 1, mismatch = -4, gapOpening = 6, gapExtension = 1, # bowtie2
) {
	# shrink anchor for small events to prevent spanning alignment
	aLength <- pmin(anchorLength, abs(gr$svLen) + 1) %na% anchorLength
	anchorSeq <- extractReferenceSequence(gr, ref, aLength, 0)
	anchorSeq <- sub(".*N", "", anchorSeq)
	# shrink anchor with Ns
	aLength <- nchar(anchorSeq)
	varseq <- extractBreakpointSequence(gr, ref, aLength)
	varseq <- sub("N.*", "", varseq)
	bpLength <- nchar(varseq) - aLength
	nonbpseq <- extractReferenceSequence(gr, ref, 0, bpLength + margin)
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


#' Converts to breakend notation
#' @param gr GRanges object.
#' @param insSeq insert sequence of the GRanges.
#' @param ref reference sequence of the GRanges.
.toVcfBreakendNotationAlt = function(gr, insSeq=gr$insSeq, ref=gr$REF) {
	assert_that(all(width(gr) == 1))
	assert_that(!is.null(insSeq))
	assert_that(all(insSeq != ""))
	assert_that(!is.null(gr$partner))
	isBreakpoint = !is.na(gr$partner)
	breakendAlt = ifelse(as.character(strand(gr)) == "+", paste0(gr$insSeq, "."), paste0(".", gr$insSeq))
	gr$partner[isBreakpoint] = names(gr)[isBreakpoint] # self partner to prevent errors
	partnergr = gr[gr$partner]
	partnerDirectionChar = ifelse(strand(partnergr) == "+", "]", "[")
	breakpointAlt = ifelse(as.character(strand(gr)) == "+",
						   paste0(ref, insSeq, partnerDirectionChar, GenomeInfoDb::seqnames(partnergr), ":", start(partnergr), partnerDirectionChar),
						   paste0(partnerDirectionChar, GenomeInfoDb::seqnames(partnergr), ":", start(partnergr), partnerDirectionChar, insSeq, ref))
	return (ifelse(isBreakpoint, breakpointAlt, breakendAlt))
}

#' Converts the given breakpoint GRanges object to VCF format in breakend
#' notation.
#'
#' @param gr breakpoint GRanges object. Can contain both breakpoint and single 
#' breakend SV records.
#' @param ... For cbind and rbind a list of VCF objects. For all other methods 
#' ... are additional arguments passed to methods. See VCF class in 
#' VariantAnnotation for more details.
#'
breakpointGRangesToVCF <- function(gr, ...) {
	if (is.null(gr$insSeq)) {
		gr$insSeq = rep("", length(gr))
	}
	nominalgr = GRanges(seqnames=GenomeInfoDb::seqnames(gr), 
	                    ranges=IRanges::IRanges(start=(end(gr) + start(gr)) / 2, 
	                                            width=1))
	if (is.null(gr$REF)) {
		gr$REF = rep("N", length(gr))
	}
	gr$ALT[is.na(gr$ALT)] = ""
	if (is.null(gr$ALT)) {
		gr$ALT = rep("", length(gr))
	}
	gr$ALT[is.na(gr$ALT)] = ""
	gr$ALT[gr$ALT == ""] = .toVcfBreakendNotationAlt(gr)[gr$ALT == ""]
	ciposstart = start(gr) - start(nominalgr)
	ciposend = end(gr) - end(nominalgr)
	vcf = VCF(rowRanges=nominalgr, collapsed=FALSE)
	fixeddf = data.frame(
		ALT=gr$ALT,
		REF=gr$REF,
		QUAL=gr$QUAL,
		FILTER=gr$FILTER)
	
	VariantAnnotation::VCF(rowRanges = GRanges(), colData = S4Vectors::DataFrame(), 
	    exptData = list(header = VCFHeader()), fixed = S4Vectors::DataFrame(), 
	    info = S4Vectors::DataFrame(), geno = S4Vectors::SimpleList(), ..., collapsed=FALSE, 
	    verbose = FALSE)

}

