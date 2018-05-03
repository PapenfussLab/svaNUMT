.dispatchPerAllele_CollapsedVCF <- function(FUN, x, singleAltOnly) {
    alt <- alt(x)
    flat <- BiocGenerics::unlist(alt, use.names=FALSE)
    res <- FUN(rep(ref(x), S4Vectors::elementNROWS(alt(x))), flat)
    lst <- relist(res, alt)
    if (singleAltOnly)
        all(lst) & S4Vectors::elementNROWS(lst) == 1
    else
        any(lst)
}
.dispatchPerAllele_ExpandedVCF <- function(FUN, x) {
    alt <- alt(x)
    flat <- BiocGenerics::unlist(alt, use.names=FALSE)
    res <- FUN(rep(ref(x), S4Vectors::elementNROWS(alt(x))), flat)
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
        grepl("]", a, fixed=TRUE) |
    	grepl(".", a, fixed=TRUE)
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
	lengthDiff <- S4Vectors::elementNROWS(ref) != IRanges::nchar(alt)
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
    a <- elementExtract(alt(vcf))
    result <- ifelse(!isStructural(vcf), 0,
		elementExtract(info(vcf)$SVLEN) %na%
		(elementExtract(info(vcf)$END) - start(rowRanges(vcf))) %na%
		(ifelse(isSymbolic(vcf), NA_integer_, IRanges::nchar(a) - IRanges::nchar(r))))
    return(result)
}

.hasSingleAllelePerRecord <- function(vcf) {
	assertthat::assert_that(is(vcf, "VCF"))
    all(S4Vectors::elementNROWS(alt(vcf)) == 1)
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
		  	.breakpointRanges(x, ...)
)

#' @param vcf VCF object
#' @param nominalPosition Determines whether to call the variant at the
#'    nominal VCF position, or to call the confidence interval (incorporating
#'    any homology present).
#' @param prefix variant name prefix to assign to unnamed variants
#' @param suffix suffix to append
#' @param info_columns VCF INFO columns to include in the GRanges object
.breakpointRanges <- function(vcf, nominalPosition=FALSE, placeholderName="svrecord", suffix="_bp", info_columns=NULL, unpartneredBreakends=FALSE) {
	vcf <- vcf[isStructural(vcf),]
	assertthat::assert_that(.hasSingleAllelePerRecord(vcf))
	# VariantAnnotation bug: SV row names are not unique
	# ensure names are defined
	if (any(duplicated(row.names(vcf)))) {
		warning("Found ", sum(duplicated(row.names(vcf))), " duplicate row names (duplicates renamed).")
	}
	if (is.null(row.names(vcf))) {
		row.names(vcf) <- paste0(placeholderName, seq_along(vcf), row.names(vcf))
	} else if (any(is.na(row.names(vcf)) | duplicated(row.names(vcf)))) {
		row.names(vcf) <- ifelse(is.na(row.names(vcf)) | duplicated(row.names(vcf)), paste0(placeholderName, seq_along(vcf)), row.names(vcf))
	}
	assertthat::assert_that(!is.null(row.names(vcf)))
	assertthat::assert_that(assertthat::noNA(row.names(vcf)))
	assertthat::assert_that(!any(duplicated(row.names(vcf))))
	gr <- rowRanges(vcf)
	gr$REF <- as.character(ref(vcf))
	gr$ALT <- as.character(elementExtract(alt(vcf), 1))
	gr$vcfId <- names(vcf)
	gr$partner <- rep(NA_character_, length(gr))
	gr$svtype <- elementExtract(info(vcf)$SVTYPE) %na%
		# hack ensure that [,2] exists even for zero record vcfs
		(stringr::str_match(c("HACK", gr$ALT), "<(.*)>")[,2][-1]) %na%
		rep(NA_character_, length(gr))
	# use the root type
	gr$svtype <- stringr::str_extract(gr$svtype, "^[^:]+")
	gr$svLen <- .svLen(vcf)
	gr$insSeq <- rep(NA_character_, length(gr))
	gr$insLen <- rep(0, length(gr))
	gr$cistartoffset <- rep(0, length(gr))
	gr$ciwidth <- rep(0, length(gr))

	for (col in info_columns) {
		mcols(gr)[[col]] <- info(vcf)[[col]]
	}
	if (!is.null(info(vcf)$HOMSEQ)) {
		seq <- elementExtract(info(vcf)$HOMSEQ, 1)
		gr$ciwidth <- ifelse(is.na(seq), gr$ciwidth, nchar(seq))
	}
	if (!is.null(info(vcf)$HOMLEN)) {
		gr$ciwidth <- elementExtract(info(vcf)$HOMLEN, 1) %na% gr$ciwidth
	}
	# have not yet factored in imprecise variant calling into ciwidth - just microhomology
	gr$HOMLEN <- gr$ciwidth

	if (!is.null(info(vcf)$CIPOS)) {
		.expectMetadataInfo(vcf, "CIPOS", 2, header.Type.Integer)
		cistartoffset <- elementExtract(info(vcf)$CIPOS, 1)
		ciendoffset <- elementExtract(info(vcf)$CIPOS, 2)
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
		if (!unpartneredBreakends) {
			commonPrefixLength <- pairwiseLCPrefix(cgr$REF, cgr$ALT, ignore.case=TRUE)
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
		}
		cgr <- NULL
		mategr <- NULL
	}
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("DEL", "INS", "DUP", "RPL")
	if (any(rows)) {
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE
		if (!unpartneredBreakends) {
			cvcf <- vcf[rows,]
			#assertthat::assert_that(!any(cgr$svtype == "DEL" & cgr$svLen > 0))
			#assertthat::assert_that(!any(cgr$svtype == "INS" & cgr$svLen < 0))
			dup <- cgr$svtype == "DUP"
			del <- cgr$svtype == "DEL"
			ins <- cgr$svtype == "INS"

			strand(cgr) <- "+"
			width(cgr) <- 1
			cgr$insLen <- ifelse(ins, abs(cgr$svLen), 0)
			if (!is.null(info(cvcf)$NTLEN)) {
				#pindel RPL
				cgr$insLen <- elementExtract(info(cvcf)$NTLEN) %na% cgr$insLen
			}
			mategr <- cgr
			strand(mategr) <- "-"
			# use end, then fall back to calculating from length
			end <- elementExtract(info(cvcf)$END, 1) %na% (start(cgr) + ifelse(ins, 0, abs(cgr$svLen)))
			if (any(is.na(end))) {
				stop(paste("Variant of undefined length: ", paste(names(cgr)[is.na(end),], collapse=", ")))
			}
			ranges(mategr) <- IRanges(start=end + ifelse(dup, 0, 1), width=1)

			cistartoffset <- elementExtract(info(cvcf)$CIEND, 1)
			ciendoffset <- elementExtract(info(cvcf)$CIEND, 2)
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
		}
		cgr <- NULL
		mategr <- NULL
	}
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("INV")
	if (any(rows)) {
		cgr1 <- gr[rows,]
		gr$processed[rows] <- TRUE
		if (!unpartneredBreakends) {
			width(cgr1) <- 1
			end <- elementExtract(info(vcf)$END[rows], 1) %na% (start(cgr1) + abs(cgr1$svLen) - 1)
			if (any(is.na(end))) {
				stop(paste("Variant of undefined length: ", paste(names(cgr1)[is.na(end),], collapse=", ")))
			}
			hasPlusBreakend <- rep(TRUE, length(cgr1))
			hasMinusBreakend <- rep(TRUE, length(cgr1))
			if (!is.null(info(vcf)$INV3)) {
				hasMinusBreakend <- !info(vcf)$INV3[rows]
			}
			if (!is.null(info(vcf)$INV5)) {
				hasPlusBreakend <- !info(vcf)$INV5[rows]
			}

			cgr2 <- cgr1
			cistartoffset <- elementExtract(info(vcf)$CIEND[rows], 1)
			ciendoffset <- elementExtract(info(vcf)$CIEND[rows], 2)
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

			outgr <- c(outgr, cgr1[hasMinusBreakend], cgr2[hasMinusBreakend], cgr3[hasPlusBreakend], cgr4[hasPlusBreakend])
		}
		cgr1 <- NULL
		cgr2 <- NULL
		cgr3 <- NULL
		cgr4 <- NULL
	}
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("BND") & str_detect(gr$ALT, "(.*)(\\[|])(.*)(\\[|])(.*)")
	if (any(rows)) {
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE
		if (!unpartneredBreakends) {
			cvcf <- vcf[rows,]

			bndMatches <- stringr::str_match(cgr$ALT, "(.*)(\\[|])(.*)(\\[|])(.*)")
			preBases <- bndMatches[,2]
			bracket <- bndMatches[,3]
			remoteLocation <- bndMatches[,4]
			postBases <- bndMatches[,6]
			strand(cgr) <- ifelse(preBases == "", "-", "+")

			cgr$partner <- NA_character_
			if (!is.null(info(cvcf)$PARID)) {
				cgr$partner <- elementExtract(info(cvcf)$PARID, 1)
			}
			if (!is.null(info(cvcf)$MATEID) & any(is.na(cgr$partner))) {
				multimates <- S4Vectors::elementNROWS(info(cvcf)$MATEID) > 1 & is.na(cgr$partner)
				cgr$partner <- ifelse(is.na(cgr$partner), elementExtract(info(cvcf)$MATEID, 1), cgr$partner)
				if (any(multimates)) {
					warning(paste("Ignoring additional mate breakends for variants", names(cgr)[multimates]))
				}
			}
			reflen <- S4Vectors::elementNROWS(cgr$REF)
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
		}
		cgr <- NULL
		mategr <- NULL
	}
	# breakend notation
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("BND")  & str_detect(gr$ALT, "^(([.].*)|(.*[.]))$")
	if (any(rows)) {
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE
		if (unpartneredBreakends) {
			cvcf <- vcf[rows,]

			strand(cgr) <- ifelse(str_sub(cgr$ALT, 1, 1) == ".", "-", "+")
			# trim anchoring base and breakend symbol
			cgr$insSeq <- str_sub(cgr$ALT, 2, str_length(cgr$ALT) - 1)
			cgr$insLen <- str_length(cgr$insSeq)
			cgr$partner <- NA_character_
			cgr$svLen <- NA_integer_
			outgr <- c(outgr, cgr)
		}
		cgr <- NULL
	}
	# TODO: Does delly write two records for a full INV?
	# DELLY TRA https://groups.google.com/forum/#!msg/delly-users/6Mq2juBraRY/BjmMrBh3GAAJ
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("TRA")
	if (any(rows)) {
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE
		if (!unpartneredBreakends) {
			cvcf <- vcf[rows,]

			if (is.null(info(cvcf)$CHR2) || any(is.na(info(cvcf)$CHR2))) {
				stop(paste("Delly variants missing CHR2:", paste(names(cgr)[is.na(info(cvcf)$CHR2)], collapse=", ")))
			}
			if (is.null(info(cvcf)$CT) || any(is.na(info(cvcf)$CT))) {
				stop(paste("Delly variants missing CT:", paste(names(cgr)[is.na(info(cvcf)$CT)], collapse=", ")))
			}
			cgr$insLen <- info(cvcf)$INSLEN %na% 0 # Delly no longer writes INSLEN to all TRA records
			width(cgr) <- 1
			mategr <- cgr
			# Hack so we can add new seqlevels if required
			seqlevels(mategr) <- unique(c(seqlevels(mategr), info(cvcf)$CHR2))
			seqnames(mategr)[seq(1, length(mategr))] <- info(cvcf)$CHR2
			ranges(mategr) <- IRanges(start=info(cvcf)$END, width=1)
			strand(cgr) <- ifelse(info(cvcf)$CT %in% c("3to3", "3to5"), "+", "-")
			strand(mategr) <- ifelse(info(cvcf)$CT %in% c("3to3", "5to3"), "+", "-")

			mcistartoffset <- elementExtract(info(cvcf)$CIEND, 1) %na% 0
			mciendoffset <- elementExtract(info(cvcf)$CIEND, 2) %na% 0
			mciwidth <- mciendoffset - mcistartoffset
			mategr$cistartoffset <- mcistartoffset
			mategr$ciwidth <- mciwidth

			names(mategr) <- paste0(names(cgr), suffix, 2)
			names(cgr) <- paste0(names(cgr), suffix, 1)
			cgr$partner <- names(mategr)
			mategr$partner <- names(cgr)
			outgr <- c(outgr, cgr, mategr)
		}
		cgr <- NULL
		mategr <- NULL
	}
	# TIGRA CTX
	rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("CTX")
	if (any(rows)) {
		# TIGRA CTX call
		cgr <- gr[rows,]
		gr$processed[rows] <- TRUE
		if (!unpartneredBreakends) {
			cvcf <- vcf[rows,]

			if (is.null(info(cvcf)$CHR2) || any(is.na(info(cvcf)$CHR2))) {
				stop(paste("TIGRA variants missing CHR2:", paste(names(cgr)[is.na(info(cvcf)$CHR2)], collapse=", ")))
			}
			width(cgr) <- 1
			mategr <- cgr
			# Hack so we can add new seqlevels if required
			seqlevels(mategr) <- unique(c(seqlevels(mategr), info(cvcf)$CHR2))
			seqnames(mategr)[seq(1, length(mategr))] <- info(cvcf)$CHR2
			ranges(mategr) <- IRanges(start=info(cvcf)$END, width=1)
			# no direction information is reported
			strand(cgr) <- "*"
			strand(mategr) <- "*"

			mcistartoffset <- elementExtract(info(cvcf)$CIEND, 1) %na% 0
			mciendoffset <- elementExtract(info(cvcf)$CIEND, 2) %na% 0
			mciwidth <- mciendoffset - mcistartoffset
			mategr$cistartoffset <- mcistartoffset
			mategr$ciwidth <- mciwidth

			names(mategr) <- paste0(names(cgr), suffix, 2)
			names(cgr) <- paste0(names(cgr), suffix, 1)
			cgr$partner <- names(mategr)
			mategr$partner <- names(cgr)
			outgr <- c(outgr, cgr, mategr)
		}
		cgr <- NULL
		mategr <- NULL
	}
	if (!all(gr$processed)) {
		stop(paste("Unrecognised format for variants", paste(names(gr)[!gr$processed], collapse=", ")))
	}
	# incorporate microhomology and confidence intervals
	if (!nominalPosition) {
		ranges(outgr) <- IRanges(start=start(outgr) + outgr$cistartoffset, width=outgr$ciwidth + 1, names=names(outgr))
	}
	outgr$processed <- NULL
	outgr$cistartoffset <- NULL
	outgr$ciwidth <- NULL
	if (!unpartneredBreakends) {
		partnerpartnerisself <- partner(outgr)$partner == names(outgr)
		if (!all(partnerpartnerisself)) {
			warning("Multiple breakends partners for a single breakend found (Ignoring all except first). StructuralVariantAnnotation does not yet support promiscuous breakpoints. Please raise a github enhancement request at https://github.com/d-cameron/StructuralVariantAnnotation/issues if this functionality is desired.")
			outgr <- outgr[partnerpartnerisself,]
		}
		# sanity check that all breakpoints partners actually exist
		haspartner <- outgr$partner %in% names(outgr)
		if (!all(haspartner)) {
			stop(paste("Sanity check failure: unpaired breakends ", paste(names(gr)[!haspartner], collapse=", ")))
		}
	} else {
		outgr$partner <- NULL
	}
	return(outgr)
}
#' Extracts unpartnered breakend structural variants as a GRanges
#'
#' See Section 5.4.9 of https://samtools.github.io/hts-specs/VCFv4.3.pdf for details.
#'
#' @export
setGeneric("breakendRanges", signature="x",
		   function(x, ...)
		   	standardGeneric("breakendRanges")
)
setMethod("breakendRanges", "VCF",
		  function(x, ...)
		  	.breakpointRanges(x, unpartneredBreakends=TRUE, ...)
)

#' @param vcf VCF object
#' @param nominalPosition Determines whether to call the variant at the
#'    nominal VCF position, or to call the confidence interval (incorporating
#'    any homology present).
#' @param prefix variant name prefix to assign to unnamed variants
#' @param suffix suffix to append
#' @param info_columns VCF INFO columns to include in the GRanges object
.breakendRanges <- function(vcf, nominalPosition=FALSE, placeholderName="svrecord", suffix="_bp", info_columns=NULL) {
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

