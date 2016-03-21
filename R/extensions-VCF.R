.dispatchPerAllele_CollapsedVCF <- function(FUN, x, singleAltOnly) {
    alt <- alt(x)
    flat <- unlist(alt, use.names=FALSE)
    res <- FUN(rep(ref(x), elementLengths(alt(x))), flat)
    lst <- relist(res, alt)
    if (singleAltOnly)
        all(lst) & elementLengths(lst) == 1
    else
        any(lst)
}
.dispatchPerAllele_ExpandedVCF <- function(FUN, x) {
    alt <- alt(x)
    flat <- unlist(alt, use.names=FALSE)
    res <- FUN(rep(ref(x), elementLengths(alt(x))), flat)
    res
}

#' Determines whether the variant is a symbolic allele
#'
#' @export
setGeneric("isSymbolic", signature="x",
           function(x, ...)
               standardGeneric("isSymbolic")
)
setMethod("isSymbolic", "CollapsedVCF",
          function(x, ..., singleAltOnly=TRUE)
              .dispatchPerAllele_CollapsedVCF(.isSymbolic, x, singleAltOnly)
)
setMethod("isSymbolic", "ExpandedVCF",
          function(x, ...)
              .dispatchPerAllele_ExpandedVCF(.isSymbolic, x)
)
.isSymbolic <- function(r, a) {
    result <- grepl("<", a, fixed=TRUE) |
        grepl("[", a, fixed=TRUE) |
        grepl("]", a, fixed=TRUE)
    return(result)
}

#' Determines whether the variant is a structural variant
#'
#' @export
setGeneric("isStructural", signature="x",
           function(x, ...)
               standardGeneric("isStructural")
)
setMethod("isStructural", "CollapsedVCF",
          function(x, ..., singleAltOnly=TRUE)
              .dispatchPerAllele_CollapsedVCF(.isStructural, x, singleAltOnly)
)
setMethod("isStructural", "ExpandedVCF",
          function(x, ...)
              .dispatchPerAllele_ExpandedVCF(.isStructural, x)
)
.isStructural <- function(ref, alt) {
	lengthDiff <- elementLengths(ref) != IRanges::nchar(alt)
	if (is(alt, "DNAStringSet")) {
		# don't break if there are no symbolic alleles in the VCF
		return(lengthDiff)
	}
    return(as.logical(
    	# exclude no-call sites
    	!is.na(alt) & alt != "<NON_REF>" &
        (lengthDiff | .isSymbolic(ref, alt))))
}


#' Returns the structural variant length of the first allele
#'
#' @param vcf VCF object
#'
.svLen <- function(vcf) {
	assertthat::assert_that(.hasSingleAllelePerRecord(vcf))
    r <- ref(vcf)
    a <- .elementExtract(alt(vcf))
    result <- ifelse(!isStructural(vcf), 0,
		.elementExtract(info(vcf)$SVLEN) %na%
		(.elementExtract(info(vcf)$END) - start(rowRanges(vcf))) %na%
		(ifelse(isSymbolic(vcf), NA_integer_, IRanges::nchar(a) - IRanges::nchar(r))))
    return(result)
}

.hasSingleAllelePerRecord <- function(vcf) {
	assertthat::assert_that(is(vcf, "VCF"))
    all(elementLengths(alt(vcf)) == 1)
}
setMethod("isStructural", "VCF",
		  function(x, ...)
		  	.dispatchPerAllele_ExpandedVCF(.isStructural, x)
)


#' Extracts the structural variants as a BreakendGRanges
#'
#' @details
#' Structural variants are converted to breakend notation.
#'
#' Due to ambiguities in the VCF specifications, structural variants
#' with multiple alt alleles are not supported.
#'
#' If HOMLEN or HOMSEQ is defined without CIPOS, it is assumed that
#' the variant position is left aligned.
#'
#' @export
setGeneric("breakpointRanges", signature="x",
		   function(x, ...)
		   		standardGeneric("breakpointRanges")
)
setMethod("breakpointRanges", "VCF",
		  function(x, ...)
		  	.breakpointRanges(x, suffix="_bp")
)
#' @param vcf VCF object
#' @param nominalPosition Determines whether to call the variant at the
#'    nominal VCF position, or to call the confidence interval (incorporating
#'    any homology present).
#' @param prefix variant name prefix to assign to unnamed variants
#' @param suffix suffix to append
.breakpointRanges <- function(vcf, nominalPosition=FALSE, placeholderName="svrecord", suffix="_bp") {
	vcf <- vcf[isStructural(vcf),]
	assertthat::assert_that(.hasSingleAllelePerRecord(vcf))
	# ensure names are defined
	if (is.null(row.names(vcf))) {
		row.names(vcf) <- paste0(placeholderName, seq_along(simple), row.names(vcf))
	} else if (any(is.na(row.names(vcf)))) {
		row.names(vcf) <- ifelse(is.na(row.names(vcf)), paste0(placeholderName, seq_along(simple), row.names(vcf)), row.names(vcf))
	}
	assertthat::assert_that(!is.null(row.names(vcf)))
	assertthat::assert_that(assertthat::noNA(row.names(vcf)))
	assertthat::assert_that(!any(duplicated(row.names(vcf))))
	gr <- rowRanges(vcf)
	gr$REF <- as.character(ref(vcf))
	gr$ALT <- as.character(.elementExtract(alt(vcf), 1))
	gr$vcfId <- names(vcf)
	gr$partner <- rep(NA_character_, length(gr))
	gr$svtype <- .elementExtract(info(vcf)$SVTYPE) %na%
		(stringr::str_match(gr$ALT, "<(.*)>")[,2]) %na%
		rep(NA_character_, length(gr))
	# use the root type
	gr$svtype <- stringr::str_extract(gr$svtype, "^[^:]+")
	gr$svLen <- .svLen(vcf)
	gr$insSeq <- rep(NA_character_, length(gr))
	gr$insLen <- rep(0, length(gr))
	gr$cistartoffset <- rep(0, length(gr))
	gr$ciwidth <- rep(0, length(gr))

	if (!is.null(info(vcf)$HOMSEQ)) {
		seq <- .elementExtract(info(vcf)$HOMSEQ, 1)
		gr$ciwidth <- ifelse(is.na(seq), gr$ciwidth, nchar(seq))
	}
	if (!is.null(info(vcf)$HOMLEN)) {
		gr$ciwidth <- .elementExtract(info(vcf)$HOMLEN, 1) %na% gr$ciwidth
	}
	if (!is.null(info(vcf)$CIPOS)) {
		.expectMetadataInfo(vcf, "CIPOS", 2, header.Type.Integer)
		cistartoffset <- .elementExtract(info(vcf)$CIPOS, 1)
		ciendoffset <- .elementExtract(info(vcf)$CIPOS, 2)
		ciwidth <- ciendoffset - cistartoffset
		gr$cistartoffset <- cistartoffset %na% gr$cistartoffset
		gr$ciwidth <- ciwidth %na% gr$ciwidth
	}
	gr$processed <- rep(FALSE, length(gr))
	outgr <- gr[FALSE,]

	rows <- !gr$processed & !isSymbolic(vcf)
	if (any(rows)) {
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE

		commonPrefixLength <- mapply(Biostrings::lcprefix, cgr$REF, cgr$ALT, USE.NAMES=FALSE)
		cgr$svLen <- nchar(cgr$ALT) - nchar(cgr$REF)
		cgr$insSeq <- subseq(cgr$ALT, start=commonPrefixLength + 1)
		cgr$insLen <- nchar(cgr$insSeq)
		start(cgr) <- start(cgr) - 1 + commonPrefixLength
		width(cgr) <- 1
		strand(cgr) <- "+"
		mategr <- cgr
		strand(mategr) <- "-"
		ranges(mategr) <- IRanges(start=start(cgr) + nchar(cgr$REF) - commonPrefixLength + 1, width=1)

		names(mategr) <- paste0(names(cgr), suffix, 2)
		names(cgr) <- paste0(names(cgr), suffix, 1)
		cgr$partner <- names(mategr)
		mategr$partner <- names(cgr)
		outgr <- c(outgr, cgr, mategr)
		cgr <- NULL
		mategr <- NULL
	}
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("DEL", "INS", "DUP")
	if (any(rows)) {
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE
		cvcf <- vcf[rows,]
		#assertthat::assert_that(!any(cgr$svtype == "DEL" & cgr$svLen > 0))
		#assertthat::assert_that(!any(cgr$svtype == "INS" & cgr$svLen < 0))
		dup <- cgr$svtype == "DUP"

		strand(cgr) <- "+"
		width(cgr) <- 1
		mategr <- cgr
		strand(mategr) <- "-"
		# use end, then fall back to calculating from length
		end <- .elementExtract(info(cvcf)$END, 1) %na% (start(cgr) + pmax(0, -cgr$svLen))
		if (any(is.na(end))) {
			stop(paste("Variant of undefined length: ", paste(names(cgr)[is.na(end),], collapse=", ")))
		}
		ranges(mategr) <- IRanges(start=end + ifelse(dup, 0, 1), width=1)
		cgr$insLen <- pmax(0, cgr$svLen)

		cistartoffset <- .elementExtract(info(cvcf)$CIEND, 1)
		ciendoffset <- .elementExtract(info(cvcf)$CIEND, 2)
		ciwidth <- ciendoffset - cistartoffset
		mategr$cistartoffset <- cistartoffset %na% mategr$cistartoffset
		mategr$ciwidth <- ciwidth %na% mategr$ciwidth


		strand(cgr)[dup] <- "-"
		strand(mategr)[dup] <- "+"

		names(mategr) <- paste0(names(cgr), suffix, 2)
		names(cgr) <- paste0(names(cgr), suffix, 1)
		cgr$partner <- names(mategr)
		mategr$partner <- names(cgr)
		outgr <- c(outgr, cgr, mategr)
		cgr <- NULL
		mategr <- NULL
	}
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("INV")
	if (any(rows)) {
		cgr1 <- gr[rows,]
		gr$processed[rows] <- TRUE

		width(cgr1) <- 1
		end <- .elementExtract(info(vcf)$END[rows], 1) %na% (start(cgr1) + abs(cgr1$svLen) - 1)
		if (any(is.na(end))) {
			stop(paste("Variant of undefined length: ", paste(names(cgr1)[is.na(end),], collapse=", ")))
		}

		cgr2 <- cgr1
		cistartoffset <- .elementExtract(info(vcf)$CIEND[rows], 1)
		ciendoffset <- .elementExtract(info(vcf)$CIEND[rows], 2)
		ciwidth <- ciendoffset - cistartoffset
		cgr2$cistartoffset <- cistartoffset %na% cgr2$cistartoffset
		cgr2$ciwidth <- ciwidth %na% cgr2$ciwidth
		cgr3 <- cgr1
		cgr4 <- cgr2

		ranges(cgr2) <- IRanges(start=end + 1, width=1)
		ranges(cgr3) <- IRanges(start=start(cgr1) - 1, width=1)
		ranges(cgr4) <- IRanges(start=end, width=1)
		strand(cgr1) <- "-"
		strand(cgr2) <- "-"
		strand(cgr3) <- "+"
		strand(cgr4) <- "+"

		names(cgr4) <- paste0(names(cgr1), suffix, 4)
		names(cgr3) <- paste0(names(cgr1), suffix, 3)
		names(cgr2) <- paste0(names(cgr1), suffix, 2)
		names(cgr1) <- paste0(names(cgr1), suffix, 1)
		cgr1$partner <- names(cgr2)
		cgr2$partner <- names(cgr1)
		cgr3$partner <- names(cgr4)
		cgr4$partner <- names(cgr3)
		outgr <- c(outgr, cgr1, cgr2, cgr3, cgr4)
		cgr1 <- NULL
		cgr2 <- NULL
		cgr3 <- NULL
		cgr4 <- NULL
	}
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("BND")
	if (any(rows)) {
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE
		cvcf <- vcf[rows,]

		bndMatches <- stringr::str_match(cgr$ALT, "(.*)(\\[|])(.*)(\\[|])(.*)")
		preBases <- bndMatches[,2]
		bracket <- bndMatches[,3]
		remoteLocation <- bndMatches[,4]
		postBases <- bndMatches[,6]
		strand(cgr) <- ifelse(preBases == "", "-", "+")

		cgr$partner <- NA_character_
		if (!is.null(info(cvcf)$PARID)) {
			cgr$partner <- .elementExtract(info(cvcf)$PARID, 1)
		}
		if (!is.null(info(cvcf)$MATEID) & any(is.na(cgr$partner))) {
			multimates <- elementLengths(info(cvcf)$MATEID) > 1 & is.na(cgr$partner)
			cgr$partner <- ifelse(is.na(cgr$partner), .elementExtract(info(cvcf)$MATEID, 1), cgr$partner)
			if (any(multimates)) {
				warning(paste("Ignoring additional mate breakends for variants", names(cgr)[multimates]))
			}
		}
		reflen <- elementLengths(cgr$REF)
		cgr$insSeq <- paste0(stringr::str_sub(preBases, reflen + 1), stringr::str_sub(postBases, end=-(reflen + 1)))
		cgr$insLen <- nchar(cgr$insSeq)

		toRemove <- is.na(cgr$partner) | !(cgr$partner %in% names(cgr))
		if (any(toRemove)) {
			warning(paste("Removing", sum(toRemove), "unpaired breakend variants", paste0(names(cgr)[toRemove], collapse=", ")))
			cgr <- cgr[!toRemove,]
		}
		mategr <- cgr[cgr$partner,]
		cgr$svLen <- ifelse(seqnames(cgr)==seqnames(mategr), abs(start(cgr) - start(mategr)) - 1, NA_integer_)
		# make deletion-like events have a -ve svLen
		cgr$svLen <- ifelse(strand(cgr) != strand(mategr) &
				((start(cgr) < start(mategr) & strand(cgr) == "+") |
				 (start(cgr) > start(mategr) & strand(cgr) == "-")),
			-cgr$svLen, cgr$svLen)
		cgr$svLen <- cgr$svLen + cgr$insLen
		outgr <- c(outgr, cgr)
		cgr <- NULL
		mategr <- NULL
	}
	# TODO: handle non-standard SVTYPE for specific callers
	# DELLY TRA
	# PINDEL RPL
	if (!all(gr$processed)) {
		stop(paste("Unrecognised format for variants", paste(names(gr)[!gr$processed,], collapse=", ")))
	}
	# incorporate microhomology and confidence intervals
	ranges(outgr) <- IRanges(start=start(outgr) + outgr$cistartoffset, width=outgr$ciwidth + 1, names=names(outgr))
	return(outgr)
}
.hasMetadataInfo <- function(vcf, field) {
	return(field %in% row.names(info(header(vcf))))
}
.expectMetadataInfo <- function(vcf, field, number, type) {
	assertthat::assert_that(.hasMetadataInfo(vcf, field))
	row <- info(header(vcf))[field,]
	assertthat::assert_that(number == row$Number)
	assertthat::assert_that(type == row$Type)
}

