#' Detecting nuclear mitochondria fusion events.
#'
#' @details
#' Nuclear mitochondrial fusion (NUMT) is a common event found in human genomes.
#' This function searches for NUMT events by identifying breakpoints supporting the fusion of
#' nuclear chromosome and mitochondrial genome. Only BND notations are supported at the current stage.
#' Possible linked nuclear insertion sites are reported by chromosome in GRanges format.
#' @param gr A GRanges object
#' @param max_ins_dist The maximum distance allowed on the reference genome between the paired insertion sites.
#' Only intra-chromosomal NUMT events are supported. Default value is 10.
#' @return A nested list of GRanges objects of candidate NUMTs.
#' @examples
#' vcf.file <- system.file("extdata", "MT.vcf", package = "NUMTDetect")
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#' gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#' numt.gr <- numtDetect(gr, max_ins_dist=20)
#' @export
#' 
numtDetect <- function(gr, max_ins_dist=10){
    assertthat::assert_that(is(gr, "GRanges"), msg = "gr should be a GRanges object")
    assertthat::assert_that(!isEmpty(gr), msg = "gr can't be empty")
    gr <- gr[seqnames(gr) %in% standardChromosomes(gr) & 
                 seqnames(partner(gr)) %in% standardChromosomes(gr)]
    pr <- breakpointgr2pairs(gr)
    numts <- pr[as.vector(xor(seqnames(S4Vectors::first(pr)) %in% c("MT", "chrM"), 
                              seqnames(S4Vectors::second(pr)) %in% c("MT", "chrM")))]
    if (length(numts)==0) {
        message("There is no NUMT event detected. Check whether 'chrM' or 'MT' is present in the VCF.")
    }else{
        #pairs objects don't ensure 
        first_MT <- seqnames(S4Vectors::first(numts)) %in% c("MT", "chrM")
        nus.gr <- S4Vectors::first(numts)
        nus.gr[first_MT] <- S4Vectors::second(numts)[first_MT]
        names(nus.gr)[as.vector(first_MT)] <- names(S4Vectors::second(numts))[as.vector(first_MT)]
        
        mts.gr <- S4Vectors::second(numts)
        mts.gr[first_MT] <- S4Vectors::first(numts)[first_MT]
        names(mts.gr)[as.vector(first_MT)] <- names(S4Vectors::first(numts))[as.vector(first_MT)]
        
        #split NU and MT grs by chromosomes
        nus <- split(nus.gr, as.vector(seqnames(nus.gr)))
        mts <- split(mts.gr, as.vector(seqnames(nus.gr)))
        #find candidate insertion sites (paired) by NU breakends close with each other (10bp)
        l <- lapply(nus, function(x) findOverlaps(x, x, maxgap = max_ins_dist, ignore.strand=TRUE))
        l <- lapply(l, function(x) x[queryHits(x)!=subjectHits(x)])
        l <- lapply(l, function(x) dplyr::as_tibble(x) %>% 
                        dplyr::mutate(label=dplyr::if_else(subjectHits>queryHits, 
                                                           paste(queryHits,subjectHits,sep = "_"), 
                                                           paste(subjectHits,queryHits,sep = "_"))) %>% 
                        dplyr::filter(!duplicated(label)) %>% dplyr::select(-label))
        if (length(unlist(l))==0){
            #no chromosome reports candidate paired insertion sites
            message("There is no NUMT event detected. Paired candidate insertion site not found.")
        }else{
            #retrieve NU breakends by results above
            n <- mapply(function(l, n) lapply(as.list(dplyr::as_tibble(t(l), .name_repair = "minimal")), 
                                              function(x) n[x]), l, nus, SIMPLIFY = FALSE)
            #retieve the corresponding MT breakends as `n`
            m <- mapply(function(l, m) lapply(as.list(dplyr::as_tibble(t(l), .name_repair = "minimal")), 
                                              function(x) m[x]), l, mts, SIMPLIFY = FALSE)
            
            ##remove empty chromosomes
            n <- n[sapply(n, function(x) length(x)>0)]
            m <- m[sapply(m, function(x) length(x)>0)] 
            
            if (sum(sapply(n, length))==0) {
                #no chromosome reports candidate paired insertion sites
                message("There is no NUMT event detected. Paired candidate insertion site not found.")
            }else{
                ##remove paired NU breakends with the same strand direction
                ##strand info is irrelevant in MT breakends, using NU's instead.
                l_strs <- lapply(n, function(x) sapply(x, function(y) as.vector(strand(y)[1]!=strand(y)[2])))
                n <- mapply(function(x, y) x[y], n, l_strs, SIMPLIFY = FALSE)
                m <- mapply(function(x, y) x[y], m, l_strs, SIMPLIFY = FALSE)
                list(NU=n, MT=m)
            }
        }
    }
}

