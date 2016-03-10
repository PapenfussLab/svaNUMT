#' Converts all symbolic structural variant calls to paired GRanges intervals.
#' 
#' @param vcf \link{VCF} object
vcfSymbolicToBreakendGRanges <- function(vcf) {
}

#' Determines whether the alternate allele is a symbolic allele
#' @export
isSymbolic <- function(vcf) {
    assert_that(is(vcf, "ExpandedVCF"))
    return(.isSymbolic(ref(vcf), alt(vcf)))
}

.isSymbolic <- function(ref, alt) {
    stringr::str_detect(alt, "[<\\[\\]]")
}
#' Determines whether the alternate allele is a structural variation
#' @export
isStructural <- function(vcf) {
    assert_that(is(vcf, "ExpandedVCF"))
    return(.isStructural(ref(vcf), alt(vcf)))
}

.isStructural <- function(ref, alt) {
    elementLengths(ref) != nchar(alt) | .isSymbolic(ref, alt)
}

svLen <- function(vcf) {
    assert_that(is(vcf, "ExpandedVCF"))
    return(.svLen(ref(vcf), alt(vcf), info(vcf)$SVTYPE, as.numeric(info(vcf)$SVLEN)))
}

.svLen <- function(ref, alt, svtype, svlen) {
    ifelse(.isStructural(ref, alt), ifelse(isSymbolic(ref, alt), svlen, elementLengths(ref) != nchar(alt)), 0)
}