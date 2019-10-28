#' Detecting nuclear mitochondria fusion events.
#'
#' @details
#' Nuclear mitochondrial fusion (NUMT) is a common event found in human genomes.
#' This function searches for NUMT events by identifying breakpoints supporting the fusion of
#' nuclear chromosome and mitochondrial genome. Only BND notations are supported at the current stage.
#' Possible linked nuclear insertion sites are reported using SV IDs in the candidatePartnerId metadata column.
#' @param gr A GRanges object
#' @param nonStandardChromosomes Whether to report insertion sites on non-standard reference 
#' chromosomes. Default value is set to FALSE.
#' @param max_ins_dist The maxium distance allowed on the reference genome between the paired insertion sites.
#' Only intra-chromosomal NUMT events are supported. Default value is 1000.
#' @return A GRanges object of possible NUMT loci.
#' @examples
#' vcf.file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#' gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#' numt.gr <- numtDetec(gr)
#' @export
numtDetec <- function(gr, nonStandardChromosomes=FALSE, max_ins_dist=1000){
    assertthat::assert_that(class(gr)=="GRanges", msg = "gr should be a GRanges object")
    assertthat::assert_that(length(gr)>0, msg = "gr can't be empty")
    if (nonStandardChromosomes==FALSE) {
        gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse", species = "Homo_sapiens")
    }
    numt.gr <- gr[!(seqnames(gr)=="chrM"|seqnames(gr)=="MT")]
    numt.gr <- numt.gr[stringr::str_match(numt.gr$ALT, "(.*)(\\[|])(.*)(:)(.+)(\\[|])(.*)")[,4] %in% c("MT", "chrM")]
    candidatePartnerId.list <- vector(mode="list", length=length(numt.gr))
    names(candidatePartnerId.list) <- names(numt.gr)
    if (length(numt.gr)>0) {
        for (i in 1:length(numt.gr)) {
            seq = seqnames(numt.gr[i])
            pos = start(numt.gr[i])
            std = strand(numt.gr[i])
            name = names(numt.gr[i])
            #candidatePartner.gr = numt.gr[seqnames(numt.gr)==seq & abs(start(numt.gr)-pos)< max_ins_dist & strand(numt.gr)!=std]
            candidatePartnerIds = names(numt.gr[seqnames(numt.gr)==seq & abs(start(numt.gr)-pos)< max_ins_dist & strand(numt.gr)!=std])
            candidatePartnerId.list[[name]] = ifelse(length(candidatePartnerIds)>0, candidatePartnerIds, NA)
            # if (length(candidatePartner.gr)>0) {
            #     candidatePartnerId = CharacterList(names(candidatePartner.gr))
            #     #candidatePartnerId = paste(candidatePartner.gr$sourceId, collapse = ",")
            #     #print(candidatePartnerId)
            #     numt.gr$candidatePartnerId[i] = candidatePartnerId
            # } else {
            #     #numt.gr$candidatePartnerId[i] = N
            #     candidatePartnerId.list[[]]
            # }
        }
        candidatePartnerId.list <- IRanges::CharacterList(candidatePartnerId.list)
        numt.gr <- c(numt.gr, gr[names(partner(gr)) %in% names(numt.gr)])
        numt.gr$candidatePartnerId <- IRanges::CharacterList()
        for (name in names(candidatePartnerId.list)) {
            numt.gr[name]$candidatePartnerId <- candidatePartnerId.list[[name]]
            numt.gr[names(partner(numt.gr))==name]$candidatePartnerId = candidatePartnerId.list[[name]]
        }
        return(numt.gr)
    }else{
        message("There is no NUMT event detected. Check whether 'chrM' or 'MT' is present in the VCF.")
    }
    #TODO: @param min_mt_len The minimum inserted mitochonrial genome length accepted. Default value is 30.
}



#' Calculating MT sequence length.
#'
#' @details
#' This function calculate the length of MT sequence length with BND notations.
#' @param bnd.start starting breakend of the MT sequence.
#' @param bnd.end ending breakend of the MT sequence.
#' @param chrM.len length of the reference MT genome.
#' @return The length of the MT sequence. When the candidate MT BNDs can't be linked as one sequence, the returned value is NA.
.mtLen <- function(bnd.start, bnd.end, chrM.len){
    bnd.start.str <- stringr::str_match(bnd.start, "(.*)(\\[|])(.*)(:)(.+)(\\[|])(.*)")
    bnd.end.str <- stringr::str_match(bnd.end, "(.*)(\\[|])(.*)(:)(.+)(\\[|])(.*)")
    assertthat::assert_that(bnd.start.str[3] %in% c("[","]"))
    assertthat::assert_that(bnd.end.str[3] %in% c("[","]"))
    assertthat::assert_that(is.numeric(bnd.start.str[6]))
    assertthat::assert_that(is.numeric(bnd.end.str[6]))
    if (bnd.start.str[3]=="[" & bnd.end.str[3]=="]") {
        if (bnd.start.str[6]<bnd.end.str[6]) {
            dist=bnd.end.str[6]-bnd.start.str[6]
        }else if (bnd.start.str[6]>=bnd.end.str[6]) {
            dist=bnd.end.str[6]-bnd.start.str[6]+chrM.len
        }
    }else if (bnd.start.str[3]=="]" & bnd.end.str[3]=="[") {
        if (bnd.start.str[6]<bnd.end.str[6]) {
            dist=chrM.len-bnd.end.str[6]+bnd.start.str[6]
        }else if (bnd.start.str[6]>=bnd.end.str[6]) {
            dist=bnd.start.str[6]-bnd.end.str[6]
        }
    }else {
        dist=NA
    }
    return(dist)
}


#' Detecting retrotranscript insertion in nuclear genomes.
#'
#' @details
#' This function searches for retroposed transcripts by identifying breakpoints supporting 
#' intronic deletions and fusions between exons and remote loci.
#' Only BND notations are supported at the current stage.
#' @param gr A GRanges object
#' @param genes TxDb object of genes. hg19 and hg38 are supported in the current version.
#' @param maxgap The maxium distance allowed on the reference genome between the paired exon boundries.
#' @param minscore The minimum proportion of intronic deletions of a transcript should be identified.
#' @return A GRangesList object, named insSite and rt, reporting breakpoints supporting insert sites and 
#' retroposed transcripts respectively. 'exon' and 'txs' in the metadata columns report exon_id and transcript_name from the 'genes' object.
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' genes <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' vcf.file <- system.file("extdata", "gridss.vcf",
#'                          package = "StructuralVariantAnnotation")
#' vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#' gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#' rt.gr <- rtDetec(gr, genes, maxgap=30, minscore=0.6)
#' @export
#' 
rtDetec <- function(gr, genes, maxgap=10, minscore=0.5){
    assertthat::assert_that(class(gr)=="GRanges", msg = "gr should be a GRanges object")
    assertthat::assert_that(length(gr)>0, msg = "gr can't be empty")
    assertthat::assert_that(class(genes)=="TxDb", msg = "genes should be a TxDb object")
    
    GenomeInfoDb::seqlevelsStyle(genes) <- GenomeInfoDb::seqlevelsStyle(gr)[1]
    genes <- GenomeInfoDb::keepSeqlevels(genes, seqlevels(genes)[1:24], pruning.mode = "coarse")
    exons <- exons(genes, columns=c("exon_id", "tx_id", "tx_name","gene_id"))
    
    hits.start <- findOverlaps(gr, exons, maxgap = maxgap, type = "start", ignore.strand = TRUE)
    hits.end <- findOverlaps(partner(gr), exons, maxgap = maxgap, type = "end", ignore.strand = TRUE)
    hits <- dplyr::inner_join(dplyr::as_tibble(hits.start), dplyr::as_tibble(hits.end), by="queryHits")
    
    if (nrow(hits)==0) {
        message("There is no NUMT event detected.")
    }else{
        txs <- mapply(intersect, exons[hits$subjectHits.x]$tx_name, exons[hits$subjectHits.y]$tx_name)
        
        rt.gr<- c(gr[hits$queryHits], partner(gr)[hits$queryHits])
        rt.gr$exon <- c(exons[hits$subjectHits.x]$exon_id, exons[hits$subjectHits.y]$exon_id)
        rt.gr$txs <- c(IRanges::CharacterList(txs), IRanges::CharacterList(txs))
        rt.gr <- rt.gr[!sapply(rt.gr$txs, rlang::is_empty)]
        #print("annotate overlapping exons")
        for (name in unique(names(rt.gr[duplicated(names(rt.gr))]))) {
            for (i in 1:length(rt.gr[names(rt.gr) == name]$txs)) {
                rt.gr[names(rt.gr) == name]$txs[[i]] <- 
                    .combineMatchingTranscripts(rt.gr, names(rt.gr[duplicated(names(rt.gr))]))[[name]]
            }
        }
        #print("merge transcript annotations in duplicated granges")
        rt.gr <- rt.gr[!duplicated(names(rt.gr))]
        #unique() and duplicated() for granges compare RANGES, not names
        rt.gr <- rt.gr[rt.gr$exon != partner(rt.gr)$exon]
        tx.rank <- .scoreByTranscripts(genes, unlist(rt.gr$txs)) 
        #dataframe of valid retro transcripts
        tx.rank <- tx.rank[tx.rank$score >= minscore,]
        
        rt.gr <- rt.gr[stringr::str_detect(unstrsplit(rt.gr$txs), paste(tx.rank$tx_name, collapse = "|"))]
        rt.gr$txs <- rt.gr$txs[mapply(stringr::str_detect, rt.gr$txs, paste(tx.rank$tx_name, collapse = "|"))]
        #select insertion site by 
        hits.start.idx <- stringr::str_detect(unstrsplit(exons[S4Vectors::subjectHits(hits.start)]$tx_name), paste(tx.rank$tx_name, collapse = "|"))
        hits.end.idx <- stringr::str_detect(unstrsplit(exons[S4Vectors::subjectHits(hits.end)]$tx_name),paste(tx.rank$tx_name, collapse = "|"))
        
        insSite.gr <- c(gr[S4Vectors::queryHits(hits.start)[hits.start.idx]], 
                        partner(gr)[S4Vectors::queryHits(hits.end)[hits.end.idx]])
        insSite.gr$exon <- c(exons[S4Vectors::subjectHits(hits.start)[hits.start.idx]]$exon_id,
                             exons[S4Vectors::subjectHits(hits.end)[hits.end.idx]]$exon_id)
        insSite.gr$txs <- c(exons[S4Vectors::subjectHits(hits.start)[hits.start.idx]]$tx_name,
                            exons[S4Vectors::subjectHits(hits.end)[hits.end.idx]]$tx_name)
        insSite.gr$txs <- insSite.gr$txs[sapply(insSite.gr$txs, function(x){x %in% tx.rank$tx_name})]
        insSite.gr <- insSite.gr[!names(insSite.gr) %in% names(rt.gr)]
        insSite.gr <- c(insSite.gr, gr[names(partner(gr)) %in% names(insSite.gr)])
        return(GRangesList(insSite = insSite.gr, rt = rt.gr))
    }
}

#' Combining matching transcripts 
#' @details
#' This is an internal function used to merge all overlapping transcripts of a breakpoint into one vector.
#' @param gr A GRanges object
#' @param names A vector of granges names.
#' @return A list of vectors. Each vector is named with the name of the corresponding granges.
.combineMatchingTranscripts <- function(gr, names){
    names <- unique(names)
    txs.list <- vector(mode="list", length=length(names))
    names(txs.list) <- names
    for (name in names) {
        #txs.list[[name]] <- name
        txs.list[[name]] <- Reduce(union, gr[names(gr) == name]$txs)
    }
    return(txs.list)
}

#' Ranking matching transcripts
#' @details
#' This is an internal function which returns overlapping transcript names with ranking scores. 
#' The ranking score is the proportion of exon-exon fusions (intronic deletion events) detected for a given transcript.
#' @param genes TxDb object of genes. hg19 and hg38 are supported in the current version.
#' @param transcripts.col A vector of transcript names.
#' @return A dataframe with two columns, tx_name and score. 
.scoreByTranscripts <- function(genes, transcripts.col){
    overlapIntron.df <- as.data.frame(table(transcripts.col)/2)
    colnames(overlapIntron.df) <- c("tx_name", "count")
    overlapIntron.df <- merge(overlapIntron.df, 
          mcols(GenomicFeatures::transcripts(genes, columns=c("tx_name","exon_rank"), 
                            filter=list(tx_name=overlapIntron.df[,1]))))
    return(data.frame(tx_name=overlapIntron.df$tx_name, 
                      score= overlapIntron.df$count / (sapply(overlapIntron.df$exon_rank, length)-1)))
}


