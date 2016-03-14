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
    gr$REF <- NULL
    gr$ALT <- NULL
    gr$vcfId <- row.names(vcf)
    gr$mateId <- row.names(gr)
    gr$svtype <- info(vcf)$SVTYPE
    gr$svLen <- svLen(vcf)
    gr$insSeq <- DNAStringSet("")
    gr$insLen <- 0
    gr$cistartoffset <- rep(0, length(gr))
    gr$ciwidth <- rep(0, length(gr))
    if (.hasMetadataInfo(vcf, "CIPOS")) {
        .expectMetadataInfo(vcf, "CIPOS", 2, header.Type.Integer)
        # Expand call position by CIPOS
        offsets <- matrix(unlist(info(vcf)$CIPOS), ncol = 2, byrow = TRUE)
        offsets[is.na(offsets)] <- 0
        gr$cistartoffset <- offsets[,1]
        gr$ciwidth <- offsets[,2] - offsets[,1]
    }
    outgr <- NULL
    
    rows <- !isSymbolic(vcf)
    if (any(rows)) {
        cgr <- gr[rows,]
        gr <- gr[!rows,]
        cvcf <- vcf[rows,]
        
        # TODO: match more than just the first base
        commonPrefixLength <- ifelse(subseq(ref(cvcf), start=1, width=1)
            == subseq(DNAStringSet(IRanges::drop(alt(cvcf))), start=1, width=1),
            1, 0)
        start(cgr) <- start(cgr) - 1 + commonPrefixLength
        strand(cgr) <- "+"
        width(cgr) <- 1
        cgr$svLen <- elementLengths(IRanges::drop(alt(cvcf))) - elementLengths(ref(cvcf))
        mategr <- cgr
        strand(mategr) <- "-"
        ranges(mategr) <- IRanges(start=start(cgr) + pmax(0, -cgr$svLen) + 1, width=1)
        cgr$insLen <- pmax(0, cgr$svLen)
        cgr$insSeq <- subseq(DNAStringSet(IRanges::drop(alt(cvcf))), start=commonPrefixLength + 1, width=cgr$insLen)
        
        cgr$mateId <- paste0(names(cgr), suffix, 2)
        names(cgr) <- paste0(names(cgr), suffix, 1)
        names(mategr) <- cgr$mateId
        mategr$mateId <- names(cgr)
        outgr <- c(outgr, cgr, mategr)
    }
    rows <- gr$svtype %in% c("DEL")
    if (any(rows)) {
        cgr <- gr[rows,]
        gr <- gr[!rows,]
        
        width(cgr) <- 1
        mategr <- cgr
        strand(cgr) <- "+"
        
        mategr <- cgr
        strand(mategr) <- "-"
        ranges(mategr) <- IRanges(start=start(cgr) + pmax(0, -cgr$svLen) + 1, width=1)
        
        cgr$mateId <- paste0(names(cgr), suffix, 2)
        names(cgr) <- paste0(names(cgr), suffix, 1)
        names(mategr) <- cgr$mateId
        mategr$mateId <- names(cgr)
        outgr <- c(outgr, cgr, mategr)
    }
    rows <- gr$svtype %in% c("INS")
    if (any(rows)) {
        cgr <- gr[rows,]
        gr <- gr[!rows,]
        
        width(cgr) <- 1
        mategr <- cgr
        strand(cgr) <- "+"
        
        mategr <- cgr
        strand(mategr) <- "-"
        ranges(mategr) <- IRanges(start=start(cgr) + pmax(0, -cgr$svLen) + 1, width=1)
        
        cgr$mateId <- paste0(names(cgr), suffix, 2)
        names(cgr) <- paste0(names(cgr), suffix, 1)
        names(mategr) <- cgr$mateId
        mategr$mateId <- names(cgr)
        outgr <- c(outgr, cgr, mategr)
    }
    rows <- gr$svtype %in% c("BND")
    if (any(rows)) {
        cgr <- gr[rows,]
        gr <- gr[!rows,]
        cvcf <- vcf[rows,]
        
        
        bndMatches <- str_match(as.character(IRanges::drop(alt(cvcf)), "(.*)(\\[|])(.*)(\\[|])(.*)")
        preBases <- bndMatches[,2]
        bracket <- bndMatches[,3]
        remoteLocation <- bndMatches[,4]
        postBases <- bndMatches[,6]
        strand(cgr) <- ifelse(preBases == "", "-", "+")
        
        
        if (!is.null(info(cvcf)$MATEID)) {
            cgr$mateId <- as.character(info(cvcf)$MATEID
        } else if (!is.null(info(vcf)$PARID)) {
            grcall[rows,]$mateIndex <- match(info(vcf)$PARID[rows], names(rowRanges(vcf)))
        }
        grcall[rows,]$untemplated <- str_length(preBases) + str_length(postBases) - str_length(as.character(rowRanges(vcf)$REF)[rows])
        if (any(rows & is.na(grcall$mateIndex))) {
            warning(paste0("Unpaired breakends ", as.character(paste(names(grcall[is.na(grcall[rows,]$mateIndex),]), collapse=", "))))
            grcall[rows & is.na(grcall$mateIndex),]$mateIndex <- seq_along(grcall)[rows & is.na(grcall$mateIndex)]
        }
        mateBnd <- grcall[grcall[rows,]$mateIndex,]
        grcall[rows,]$size <- ifelse(seqnames(mateBnd)==seqnames(grcall[rows,]), abs(start(grcall[rows,]) - start(mateBnd)) - 1 + grcall[rows,]$untemplated, NA_integer_)
    }
    if (length(gr) > 0) {
       stop(paste("Unrecognised format for variants ", paste(names(gr), collapse=", ")))
    }
    return(outgr)
}

mateBreakend <- function(begr) {
    assert_that(begr$mateId %in% row.names(begr))
    return(begr[begr$mateId,])
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

