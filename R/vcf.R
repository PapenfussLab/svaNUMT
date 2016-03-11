# VariantAnnotation extensions

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
    as.logical(alt != "<NON_REF>" & #gVCF no-call site
        (elementLengths(ref) != IRanges::nchar(alt) | .isSymbolic(ref, alt)))
}


#' Returns the structural variant length
#' 
#' @param vcf VCF object
#' 
#' @export
svLen <- function(vcf) {
    assert_that(.hasSingleAllelePerRecord(vcf))
    r <- ref(vcf)
    a <- alt(vcf)
    svlen <- as.integer(info(vcf)$SVLEN)
    return(ifelse(!isStructural(vcf), 0,
        ifelse(isSymbolic(vcf), svlen, as.integer(IRanges::nchar(a)) - elementLengths(r))))
}

.hasSingleAllelePerRecord <- function(vcf) {
    all(elementLengths(alt(vcf)) == 1)
}

#' Converts all structural variant calls to paired GRanges intervals.
#' 
#' @param vcf VCF object
#' 
#' @export
vcfSVToBreakendGRanges <- function(vcf) {
    assert_that(is(vcf, "ExpandedVCF"))
    vcf <- vcf[isStructural(vcf),]
}
