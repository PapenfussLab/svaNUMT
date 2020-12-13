#' Converting breakpoint GRanges to BEDPE-like dataframe
#' @details
#' \code{breakpointgr2bedpe} converts a breakpoint GRanges to a BEDPE-formatted
#' dataframe. The BEDPE format consists of two sets of genomic loci, optional
#' columns of name, score, strand1, strand2 and any user-defined fields.
#' See \url{https://bedtools.readthedocs.io/en/latest/content/general-usage.html}
#' for more details on the BEDPE format.
#' @param gr A GRanges object.
#' @return A BEDPE-formatted data frame.
#' @examples
#' #coverting a GRanges object to BEDPE-like dataframe
#' vcf.file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#' gr <- breakpointRanges(vcf)
#' breakpointgr2bedpe(gr)
#' @export
breakpointgr2bedpe <- function(gr){
	bedpe <- data.frame(
		chrom1=GenomeInfoDb::seqnames(gr),
		start1=start(gr) - 1,
		end1=end(gr),
		chrom2=GenomeInfoDb::seqnames(partner(gr)),
		start2=start(partner(gr)) - 1,
		end2=end(partner(gr)),
		name=names(gr),
		partner.name=names(partner(gr)),
		score=gr$QUAL,
		strand1=strand(gr),
		strand2=strand(partner(gr)))
	bedpe <- bedpe[(as.character(bedpe$chrom1)<as.character(bedpe$chrom2)) |
		  	(bedpe$chrom1==bedpe$chrom2 & bedpe$start1<bedpe$start2),-c(8)]
	return(bedpe)
}
