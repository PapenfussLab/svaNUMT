#' Detecting nuclear mitochondria fusion events.
#'
#' @details
#' Nuclear mitochondrial fusion (NUMT) is a common event found in human genomes.
#' This function searches for NUMT events by identifying breakpoints supporting the fusion of
#' nuclear chromosome and mitochondrial genome. Only BND notations are supported at the current stage.
#' Possible linked nuclear insertion sites are reported by chromosome in GRanges format.
#' @param gr A GRanges object
#' @inheritParams numtDetect_MT
#' @inheritParams numtDetect_insseq
#' @inheritParams numtDetect_known
#' @return A nested list of GRanges objects of candidate NUMTs.
#' @examples
#' vcf.file <- system.file("extdata", "MT.vcf", package = "svaNUMT")
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#' gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#' numtS <- readr::read_table(system.file("extdata", "numtS.txt", package = "svaNUMT"), col_names = FALSE)
#' colnames(numtS) <- c("bin", "seqnames", "start", "end", "name", "score", "strand")
#' numtS <- `seqlevelsStyle<-`(GRanges(numtS), "NCBI")
#' genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' genomeMT <- genome$chrMT
#' numt.gr <- numtDetect(gr, numtS, genomeMT, max_ins_dist=20)
#' @export
#' 
#' 
numtDetect <- function(gr, numtS, genomeMT, max_ins_dist=10, maxgap_numtS=10, 
                       min_len=20, min.Align=0.8){
    list(MT=numtDetect_MT(gr, max_ins_dist = max_ins_dist),
         known=numtDetect_known(gr, numtS, max_ins_dist=max_ins_dist, maxgap_numtS=maxgap_numtS),
         insSeq=numtDetect_insseq(gr, genomeMT, min_len=min_len, min.Align=min.Align))
}


