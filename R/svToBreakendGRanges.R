#' Converts all structural variant calls to paired GRanges intervals.
#' 
#' @details
#' All rows must be named
#' 
#' If HOMLEN or HOMSEQ is defined without CIPOS, it is assumed that
#' the variant position is left aligned
#' 
#' @param vcf VCF object
#' @param prefix variant name prefix to assign to unnamed variants
#' @param suffix suffix to append 
#' 
#' @export
svToBreakendGRanges <- function(vcf, suffix="_bp") {
    assert_that(.hasSingleAllelePerRecord(vcf))
    vcf <- vcf[isStructural(vcf),]
    # ensure names are defined
    if (is.null(row.names(vcf))) {
        row.names(vcf) <- paste0("svrecord", seq_along(simple), row.names(vcf))
    } else if (any(is.na(row.names(vcf)))) {
        row.names(vcf) <- ifelse(is.na(row.names(vcf)), paste0("svrecord", seq_along(simple), row.names(vcf)), row.names(vcf))
    }
    assert_that(!is.null(row.names(vcf)))
    assert_that(noNA(row.names(vcf)))
    assert_that(!any(duplicated(row.names(vcf))))
    gr <- rowRanges(vcf)
    gr$REF <- as.character(ref(vcf))
    gr$ALT <- as.character(.elementExtract(alt(vcf), 1))
    gr$vcfId <- row.names(vcf)
    gr$partner <- rep(NA_character_, length(gr))
    gr$svtype <- rep(NA_integer_, length(gr))
    if (!is.null(info(vcf)$SVTYPE)) {
        gr$svtype <- info(vcf)$SVTYPE
    }
    gr$svLen <- svLen(vcf)
    gr$insSeq <- rep(NA_character_, length(gr))
    gr$insLen <- rep(0, length(gr))
    gr$cistartoffset <- rep(0, length(gr))
    gr$ciwidth <- rep(0, length(gr))
    if (!is.null(info(vcf)$CIPOS)) {
        .expectMetadataInfo(vcf, "CIPOS", 2, header.Type.Integer)
        gr$cistartoffset <- .elementExtract(info(vcf)$CIPOS, 1)
        gr$ciwidth <- .elementExtract(info(vcf)$CIPOS, 2) -
            .elementExtract(info(vcf)$CIPOS, 1)
    }
    gr$processed <- rep(FALSE, length(gr))
    outgr <- gr[FALSE,]
    
    rows <- !gr$processed & !isSymbolic(vcf)
    if (any(rows)) {
        cgr <- gr[rows,]
        gr$processed[rows] <- TRUE
        
        commonPrefixLength <- mapply(Biostrings::lcprefix, cgr$REF, cgr$ALT, USE.NAMES=FALSE)
        cgr$svLen <- nchar(cgr$ALT) - nchar(cgr$REF)
        cgr$insSeq <- subseq(cgr$ALT, start=commonPrefixLength + 1, width=cgr$insLen)
        cgr$insLen <- nchar(cgr$insSeq)
        start(cgr) <- start(cgr) - 1 + commonPrefixLength
        width(cgr) <- 1
        strand(cgr) <- "+"
        mategr <- cgr
        strand(mategr) <- "-"
        ranges(mategr) <- IRanges(start=start(cgr) + 1 + pmax(-cgr$svLen, 0), width=1)
                
        cgr$partner <- paste0(names(cgr), suffix, 2)
        names(cgr) <- paste0(names(cgr), suffix, 1)
        names(mategr) <- cgr$partner
        mategr$partner <- names(cgr)
        outgr <- c(outgr, cgr, mategr)
    }
    rows <- !gr$processed & !is.na(gr$svtype) & gr$svtype %in% c("DEL", "INS")
    if (any(rows)) {
        cgr <- gr[rows,]
        gr$processed[rows] <- TRUE
        assert_that(!any(cgr$svtype == "DEL" & cgr$svLen > 0))
        assert_that(!any(cgr$svtype == "INS" & cgr$svLen < 0))
        
        strand(cgr) <- "+"
        width(cgr) <- 1
        mategr <- cgr
        strand(mategr) <- "-"
        ranges(mategr) <- IRanges(start=start(cgr) + pmax(0, -cgr$svLen) + 1, width=1)
    	cgr$insLen <- pmax(0, cgr$svLen)
        
        cgr$partner <- paste0(names(cgr), suffix, 2)
        names(cgr) <- paste0(names(cgr), suffix, 1)
        names(mategr) <- cgr$partner
        mategr$partner <- names(cgr)
        outgr <- c(outgr, cgr, mategr)
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
            cgr$partner <- ifelse(is.na(cgr$partner), .elementExtract(info(cvcf)$MATEID, 1), cgr$partner)
            multimates <- elementLengths(info(cvcf)$partner) > 1 & is.na(cgr$partner)
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
        cgr$svLen <- ifelse(seqnames(cgr)==seqnames(mategr), abs(start(cgr) - start(mategr)) - 1 + cgr$insLen, NA_integer_)
        outgr <- c(outgr, cgr)
    }
    
    if (!all(gr$processed)) {
       stop(paste("Unrecognised format for variants", paste(names(gr)[!gr$processed,], collapse=", ")))
    }
    return(outgr)
}

partner <- function(begr) {
    assert_that(all(begr$partner %in% names(begr)))
    return(begr[begr$partner,])
}
.hasMetadataInfo <- function(vcf, field) {
    return(field %in% row.names(info(header(vcf))))
}
.expectMetadataInfo <- function(vcf, field, number, type) {
    assert_that(.hasMetadataInfo(vcf, field))
    row <- info(header(vcf))[field,]
    assert_that(number == row$Number)
    assert_that(type == row$Type)
}

