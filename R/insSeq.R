#' Detecting nuclear mitochondria fusion events from unmapped insertion sequences.
#'
#' @details
#' This function looks for NUMTs which the insertion MT sequences come from 
#' insertion sequences reported by SV callers.
#' @param gr A GRanges object
#' @param genomeMT A genome object of the mitochondria.
#' @param min_len The minimum length allowed of the insertion sequences.
#' Default value is 20.
#' @param min.Align The minimum alignment score allowed between the insertion 
#' sequence and MT genome.
#' @keywords internal
#' @return A nested list of GRanges objects of candidate NUMTs.
numtDetect_insseq <- function(gr, genomeMT, min_len=20, min.Align=0.8){
    assertthat::assert_that(is(gr, "GRanges"), msg = "gr should be a GRanges object")
    assertthat::assert_that(is(numtS, "GRanges"), msg = "numtS should be a GRanges object")
    assertthat::assert_that(!isEmpty(gr), msg = "gr can't be empty")
    gr <- gr[seqnames(gr) %in% standardChromosomes(gr) & 
                 seqnames(partner(gr)) %in% standardChromosomes(gr)]
    pr <- breakpointgr2pairs(gr)
    numts <- pr[as.vector(nchar(S4Vectors::first(pr)$insSeq)>=min_len)]
    numts <- numts[as.vector(!(seqnames(S4Vectors::first(numts)) %in% c("MT", "chrM") | 
                                   seqnames(S4Vectors::second(numts)) %in% c("MT", "chrM")))]
    #isEmpty() is not defined for objects of class Pairs
    if (length(numts)==0) {
        message("There is no NUMT event detected by insertion sequences.")
    }else{
        insSeqs <- S4Vectors::first(numts)$insSeq
        scores <- vapply(insSeqs, function(x) seqAlignment.score(seq=x, genome = genomeMT), 
                         FUN.VALUE = numeric(1), USE.NAMES = FALSE)
        i <- as.vector(scores >= min.Align)
        if (sum(i)==0) {
            message("There is no NUMT event detected satisfying the min.Align score.")
        }else{
            numt.insseq <- c(S4Vectors::first(numts)[i], 
                             S4Vectors::second(numts)[i])
            numt.insseq$name <- rep(S4Vectors::first(numts)[i]$sourceId, 2)
            numt.insseq$score <- rep(scores[i], 2)
            numt.insseq <- split(numt.insseq, as.vector(seqnames(numt.insseq)))
            lapply(numt.insseq, function(x) split(x, as.vector(x$name)))
            
        }
    }
}

#' Calculating the alignment score between a DNA sequence and target genome.
#'
#' @details
#' This function calculates the alignment score between a DNA sequence and target genome.
#' @param seq A string of DNA sequence.
#' @param genome An XString of the target genome.
#' @keywords internal
#' @return A alignment score between a DNA sequence and target genome.
seqAlignment.score <- function(seq, genome){
    mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1)
    seq <- DNAString(seq) #NA seq returns a 0-letter object
    s1 <- pairwiseAlignment(pattern = seq, subject = genome, 
                            scoreOnly=TRUE,  type="local", substitutionMatrix = mat,
                            gapOpening=0, gapExtension=1)/length(seq)
    s2 <- pairwiseAlignment(pattern = seq, subject = reverseComplement(genome), 
                            scoreOnly=TRUE,  type="local", substitutionMatrix = mat,
                            gapOpening=0, gapExtension=1)/length(seq)
    max(s1, s2)
}
