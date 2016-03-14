#' testthat helper utility to locate files used
#' for package tests
.testfile <- function(filename, location="extdata") {
    if (file.exists(filename)) return(filename)
    f <- system.file(location, filename, package="StructuralVariantAnnotation")
    if (!file.exists(f)) {
        f <- file.path(getwd(), "inst", location, filename)
    }
    assert_that(file.exists(f))
    return(f)
}
#' loads a VCF containing the given records
#' @param record string vector of record to write
.testrecord <- function(record) {
    filename=tempfile(fileext=".vcf")
    write(paste0(c(
        "##fileformat=VCFv4.2",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=INV,Description=\"Inversion\">",
        "##ALT=<ID=TRA,Description=\"Translocation\">",
        "##ALT=<ID=INS,Description=\"Insertion\">",
        "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">",
        "##contig=<ID=chr1,length=249250621>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        record), collapse="\n"),
        file=filename)
    vcf <- readVcf(.testfile(filename), "")
    file.remove(filename)
    return(vcf)
}