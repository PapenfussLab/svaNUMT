#' testthat helper utility to locate files used
#' for package tests
.testfile <- function(filename, location="extdata") {
    if (file.exists(filename)) return(filename)
    f <- system.file(location, filename, package="StructuralVariantAnnotation")
    if (!file.exists(f)) {
        f <- file.path(getwd(), "inst", location, filename)
    }
    assertthat::assert_that(file.exists(f))
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
        "##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">",
        "##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">",
        "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">",
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=PARID,Number=1,Type=String,Description=\"ID of partner breakend\">",
		"##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">",
    	"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    	"##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">",
    	"##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">",
    	"##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">",
        "##contig=<ID=chr1,length=249250621>",
    	"##contig=<ID=chr12>",
    	"##contig=<ID=chrM,length=16571>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        record), collapse="\n"),
        file=filename)
    vcf <- readVcf(.testfile(filename), "")
	# readVcf holds a file handle open so we won't be able to delete the
	# file even if we tried
    #file.remove(filename)
    return(vcf)
}
