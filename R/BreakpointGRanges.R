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
#' @param selfPartnerSingleBreakends treat single breakends as their own partner.
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
partner <- function(gr, selfPartnerSingleBreakends=FALSE) {
  .assertValidBreakpointGRanges(gr, allowSingleBreakends=selfPartnerSingleBreakends)
  return(gr[ifelse(selfPartnerSingleBreakends & is.na(gr$partner), names(gr), gr$partner),])
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
  .assertValidBreakpointGRanges(query)
  .assertValidBreakpointGRanges(subject)
  pquery = partner(query)
  squery = partner(subject)
  localhits = findOverlaps(query, subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand)
  remotehits = findOverlaps(pquery, squery, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=ignore.strand)
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
    callwidth <- .distance(query, pquery)
    truthwidth <- .distance(subject, squery)
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
    remotebperror <- .distance(pquery[S4Vectors::queryHits(hits)], squery[S4Vectors::subjectHits(hits)])$min
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
#' The best breakpoint is considered to be the one with the highest QUAL score.
#' @param breakpointScoreColumn Query column defining a score for determining which query breakpoint
#' is considered the best when countOnlyBest=TRUE.
#' @examples
#' truth_vcf = VariantAnnotation::readVcf(system.file("extdata", "na12878_chr22_Sudmunt2015.vcf", package = "StructuralVariantAnnotation"))
#' crest_vcf = VariantAnnotation::readVcf(system.file("extdata", "na12878_chr22_crest.vcf", package = "StructuralVariantAnnotation"))
#' caller_bpgr = breakpointRanges(crest_vcf)
#' caller_bpgr$true_positive = countBreakpointOverlaps(caller_bpgr, breakpointRanges(truth_vcf),
#'   maxgap=100, sizemargin=0.25, restrictMarginToSizeMultiple=0.5, countOnlyBest=TRUE)
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
      dplyr::arrange(dplyr::desc(.data$QUAL), .data$queryHits) %>%
      dplyr::distinct(.data$subjectHits, .keep_all=TRUE) %>%
      dplyr::group_by(.data$queryHits) %>%
      dplyr::summarise(n=dplyr::n())
  }
  hitscounts[hits$queryHits] <- hits$n
  return(hitscounts)
}

#' Converts a breakpoint GRanges object to a Pairs object
#' @param bpgr breakpoint GRanges object
#' @param writeQualAsScore write the breakpoint GRanges QUAL field as the score
#' fields for compatibility with BEDPE rtracklayer export
#' @param writeName write the breakpoint GRanges QUAL field as the score
#' fields for compatibility with BEDPE rtracklayer export
#' @param bedpeName function that returns the name to use for the breakpoint.
#' Defaults to the sourceId, name column, or row names (in that priority) of
#' the first breakend of each pair.
#' @param firstInPair function that returns TRUE for breakends that are considered
#' the first in the pair, and FALSE for the second in pair breakend. By default,
#' the first in the pair is the breakend with the lower ordinal in the breakpoint
#' GRanges object.
#' @examples
#' vcf.file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
#' bpgr <- breakpointRanges(VariantAnnotation::readVcf(vcf.file))
#' pairgr <- breakpointgr2pairs(bpgr)
#' rtracklayer::export(pairgr, con="example.bedpe")
#' @return Pairs GRanges object suitable for export to BEDPE by rtracklayer
#' @rdname pairs2breakpointgr
#' @export
breakpointgr2pairs <- function(
  bpgr,
  writeQualAsScore=TRUE,
  writeName=TRUE,
  bedpeName = NULL,
  firstInPair = NULL) {
  .assertValidBreakpointGRanges(bpgr, "Cannot convert breakpoint GRanges to Pairs: ", allowSingleBreakends=FALSE)
  
  if (is.null(bedpeName)) {
    bedpeName = function(gr) { (gr$sourceId %null% gr$name) %null% names(gr) }
  }
  if (is.null(firstInPair)) {
    firstInPair = function(gr) { seq_along(gr) < match(gr$partner, names(gr)) }
  }
  isFirst = firstInPair(bpgr)
  pairgr = S4Vectors::Pairs(bpgr[isFirst], partner(bpgr)[isFirst])
  if (writeName) {
    S4Vectors::mcols(pairgr)$name = bedpeName(S4Vectors::first(pairgr))
  }
  if (writeQualAsScore) {
    S4Vectors::mcols(pairgr)$score = S4Vectors::first(pairgr)$QUAL
  }
  return(pairgr)
}
.assertValidBreakpointGRanges <- function(bpgr, friendlyErrorMessage="", allowSingleBreakends=TRUE) {
  if (is.null(names(bpgr))) {
    stop(paste0(friendlyErrorMessage, "Breakpoint GRanges require names"))
  }
  if (any(is.na(names(bpgr)))) {
    stop(paste0(friendlyErrorMessage, "Breakpoint GRanges names cannot be NA"))
  }
  if (any(duplicated(names(bpgr)))) {
    stop(paste0(friendlyErrorMessage, "Breakpoint GRanges names cannot duplicated"))
  }
  if (!allowSingleBreakends & any(is.na(bpgr$partner))) {
    stop(paste0(friendlyErrorMessage, "Breakpoint GRanges contains single breakends"))
  }
  if (any(duplicated(bpgr$partner) & !is.na(bpgr$partner))) {
    stop(paste0(friendlyErrorMessage,
                "Multiple breakends with the sample partner identified. ",
                "Breakends with multiple partners not currently supported by Breakpoint GRanges."))
  }
  else if (!all(is.na(bpgr$partner) | (bpgr$partner %in% names(bpgr) & names(bpgr) %in% bpgr$partner))) {
    stop(paste0(friendlyErrorMessage,
                "Unpartnered breakpoint found. ",
                "All breakpoints must contain a partner in the breakpoint GRanges."))
  }
}
#' Converts a BEDPE Pairs containing pairs of GRanges loaded using to a breakpoint GRanges object.
#' @details
#' Breakpoint-level column names will override breakend-level column names.
#' @param pairs a Pairs object consisting of two parallel genomic loci.
#' @param placeholderName prefix to use to ensure each entry has a unique ID.
#' @param firstSuffix first in pair name suffix to ensure breakend name uniqueness
#' @param secondSuffix second in pair name suffix to ensure breakend name uniqueness
#' @param nameField Fallback field for row names if the Pairs object does not contain any names.
#' BEDPE files loaded using rtracklayer use the "name" field.
#' @param renameScoreToQUAL renames the 'score' column to 'QUAL'.
#' Performing this rename results in a consistent variant quality score column
#' name for variant loaded from BEDPE and VCF.
#' @examples
#' bedpe.file <- system.file("extdata", "gridss.bedpe", package = "StructuralVariantAnnotation")
#' bedpe.pairs <- rtracklayer::import(bedpe.file)
#' bedpe.bpgr <- pairs2breakpointgr(bedpe.pairs)
#' @return Breakpoint GRanges object.
#' @export
pairs2breakpointgr <- function(
		pairs,
		placeholderName="bedpe",
		firstSuffix="_1", secondSuffix="_2",
		nameField="name",
		renameScoreToQUAL=TRUE) {
	n <- names(pairs)
	if (is.null(n)) {
		# BEDPE uses the "name" field
		if (nameField %in% names(S4Vectors::mcols(pairs))) {
			n <- S4Vectors::mcols(pairs)[[nameField]]
			mcols(pairs)$sourceId <- n
		} else {
			n <- rep(NA_character_, length(pairs))
		}
	}
	# ensure row names are unique
	n <- ifelse(is.na(n) | n == "" | n =="." | duplicated(n), paste0(placeholderName, seq_along(n)), n)
	#
	gr <- c(S4Vectors::first(pairs), S4Vectors::second(pairs))
	names(gr) <- c(paste0(n, firstSuffix), paste0(n, secondSuffix))
	gr$partner <- c(paste0(n, secondSuffix), paste0(n, firstSuffix))
	for (col in names(S4Vectors::mcols(pairs))) {
		if (col %in% nameField) {
			# drop columns we have processed
		} else {
			S4Vectors::mcols(gr)[[col]] <- S4Vectors::mcols(pairs)[[col]]
		}
	}
	if (renameScoreToQUAL) {
		names(mcols(gr))[which(names(mcols(gr)) == "score")] <- "QUAL"
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
#' @return Breakpoint sequence around the variant position.
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
#' @return Reference sequence around the breakpoint position.
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
#' @return A constricted GRanges object.
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
#' @return A dataframe containing the length of inexact homology between the 
#' breakpoint sequence and the reference.
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
#' @return breakendAlt or breakpointAlt depending on whether the variant is partnered.
.toVcfBreakendNotationAlt = function(gr, insSeq=gr$insSeq, ref=gr$REF) {
	assertthat::assert_that(all(width(gr) == 1))
	assertthat::assert_that(!is.null(insSeq))
	assertthat::assert_that(all(insSeq != ""))
	assertthat::assert_that(!is.null(gr$partner))
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
#' @return A VCF object.
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
#' Type of simplest explaination of event. Possible types are:
#' | Type | Description |
#' | BND | Single breakend |
#' | CTX | Interchromosomal translocation |
#' | INV | Inversion. Note that both ++ and -- breakpoint will be classified as inversion regardless of whether the matching breakpoint actually exists |
#' | DUP | Tandem duplication |
#' | INS | Insertion |
#' | DEL | Deletion |
#' 
#' @param gr breakpoint GRanges object
#' @param insertionLengthThreshold portion of inserted bases compared to total event size to be classified as an insertion. For example, a 5bp deletion with 5 inserted bases will be classified as an INS event.
#' @return Type of simplest explaination of event
#' @export
simpleEventType <- function(gr, insertionLengthThreshold=0.5) {
	if (is.null(gr$partner)) {
		gr$partner = rep(NA_character_, length(gr))
	}
	pgr = partner(gr, selfPartnerSingleBreakends=TRUE)
	return(
		ifelse(is.na(gr$partner), "BND", 
			ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
				ifelse(strand(gr) == strand(pgr), "INV",
					ifelse(gr$insLen >= abs(simpleEventLength(gr)) * insertionLengthThreshold, "INS", # TODO: improve classification of complex events
						ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
							"DUP"))))))
}
#' Length of event if interpreted as an isolated breakpoint.
#' @param gr breakpoint GRanges object
#' @return Length of the simplest explaination of this breakpoint/breakend.
#' @export
simpleEventLength <- function(gr) {
	if (is.null(gr$partner)) {
		gr$partner = rep(NA_character_, length(gr))
	}
	pgr = partner(gr, selfPartnerSingleBreakends=TRUE)
	return(
		ifelse(seqnames(gr) != seqnames(pgr) | as.logical(strand(gr) == strand(pgr) | is.na(gr$partner)), NA_integer_,
			gr$insLen + 1 + ifelse(as.logical(strand(gr) == "+"), start(gr) - start(pgr), start(pgr) - start(gr))))
}
#' Finds duplication events that are reported as inserts.
#' As sequence alignment algorithms do no allow backtracking, long read-based
#' variant callers will frequently report small duplication as insertion events.
#' Whilst both the duplication and insertion representations result in the same
#' sequence, this representational difference is problematic when comparing
#' variant call sets.
#' 
#' WARNING: this method does not yet check that the inserted sequence actually matched the duplicated sequence.
#' @param query a breakpoint GRanges object
#' @param subject a breakpoint GRanges object
#' @param maxgap maximum distance between the insertion position and the duplication
#' @param maxsizedifference maximum size difference between the duplication and insertion.
#' @return Hits object containing the ordinals of the matching breakends
#' in the query and subject 
#' @export
findInsDupOverlaps <- function(query, subject, maxgap=-1L, maxsizedifference=0L) {
	.assertValidBreakpointGRanges(query)
	.assertValidBreakpointGRanges(subject)
	query$ordinal = seq_len(length(query))
	subject$ordinal = seq_len(length(subject))
	query$set = simpleEventType(query)
	query$sel = simpleEventLength(query)
	subject$set = simpleEventType(subject)
	subject$sel = simpleEventLength(subject)
	pquery = partner(query)
	psubject = partner(subject)
	query$isLowBreakend = start(query) < start(pquery) | (start(query) == start(pquery) & query$ordinal < pquery$ordinal)
	subject$isLowBreakend = start(subject) < start(psubject) | (start(subject) == start(psubject) & subject$ordinal < psubject$ordinal)
	
	qins_to_sdup = .findOverlaps_queryIns_subjectDup(query, subject, psubject, maxgap=maxgap, maxsizedifference=maxsizedifference)
	sins_to_qdup = .findOverlaps_queryIns_subjectDup(subject, query, pquery, maxgap=maxgap, maxsizedifference=maxsizedifference)
	lowhits = data.frame(
		qhits=c(qins_to_sdup$queryHits, sins_to_qdup$subjectHits),
		shits=c(qins_to_sdup$subjectHits, sins_to_qdup$queryHits))
	# add upper to upper match
	bothhits = Hits(
		from=c(lowhits$qhits, pquery$ordinal[lowhits$qhits]),
		to=c(lowhits$shits, psubject$ordinal[lowhits$shits]),
		nLnode=length(query),
		nRnode=length(subject))
	return(bothhits)
}
.findOverlaps_queryIns_subjectDup <- function(query, subject, psubject , maxgap=-1L, maxsizedifference=0L) {
	subject$HighEndPosition = end(psubject)
	subject = subject[subject$set == "DUP" & subject$isLowBreakend]
	end(subject) = subject$HighEndPosition
	# expand by one since insertion can preceed, succeed, or be in the middle of the dup
	start(subject) = start(subject) - 1
	query = query[query$set == "INS" & query$isLowBreakend]
	hits = findOverlaps(query, subject, maxgap=maxgap, ignore.strand=TRUE)
	hits = hits[abs(query$sel[queryHits(hits)] - subject$sel[subjectHits(hits)]) <= maxsizedifference]
	# TODO: filter by 
	# Translate back to ordinals of what was passed in to us
	return(data.frame(
		queryHits=query$ordinal[queryHits(hits)],
		subjectHits=subject$ordinal[subjectHits(hits)]))
}
