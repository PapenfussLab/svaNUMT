rtDetec <- function(gr, genes, maxgap=100, minscore=0.3){
    #message("26.01.2020")
    #check args
    assertthat::assert_that(class(gr)=="GRanges", msg = "gr should be a GRanges object")
    assertthat::assert_that(length(gr)>0, msg = "gr can't be empty")
    assertthat::assert_that(class(genes)=="TxDb", msg = "genes should be a TxDb object")
    
    #prepare annotation exons
    GenomeInfoDb::seqlevelsStyle(genes) <- GenomeInfoDb::seqlevelsStyle(gr)[1]
    genes <- GenomeInfoDb::keepSeqlevels(genes, seqlevels(genes)[1:24], pruning.mode = "coarse")
    exons <- exons(genes, columns=c("exon_id", "tx_id", "tx_name","gene_id"))
    
    #------------------------
    #testing only
    gr <- breakpointRanges(manta)
    #------------------------
    
    #find exon-SV overlaps:
    hits.start <- findOverlaps(gr, exons, maxgap = maxgap, type = "start", ignore.strand = TRUE)
    hits.end <- findOverlaps(partner(gr), exons, maxgap = maxgap, type = "end", ignore.strand = TRUE)
    
    # 1.return breakpoints overlapping with exons on both ends (>=2 exons)
    hits <- dplyr::inner_join(dplyr::as_tibble(hits.start), dplyr::as_tibble(hits.end), by="queryHits")
    #mcols(exons)[hits$subjectHits.x, "gene_id"] == mcols(exons)[hits$subjectHits.y, "gene_id"]
    same.tx <- sapply(Reduce(intersect, list(mcols(exons)[hits$subjectHits.x, 'tx_id'], 
                                             mcols(exons)[hits$subjectHits.y, 'tx_id'])),length)!=0
    hits.tx <- hits[same.tx,]
    
    # 2.return breakpoints of insertionSite-exon 
    hits.insSite <- hits[!same.tx,] %>% 
        bind_rows(.,anti_join(dplyr::as_tibble(hits.start), dplyr::as_tibble(hits.end), by='queryHits')) %>% 
        bind_rows(., anti_join(dplyr::as_tibble(hits.end), dplyr::as_tibble(hits.start), by='queryHits'))
            
        
    
    if (nrow(hits.tx)+nrow(hits.insSite)==0) {
        message("There is no retroposed gene detected.")
        return(GRanges())
    }else{
        txs <- mapply(intersect, exons[hits.tx$subjectHits.x]$tx_name, exons[hits.tx$subjectHits.y]$tx_name)
        
        rt.gr<- c(gr[hits.tx$queryHits], partner(gr)[hits.tx$queryHits])
        rt.gr$exon <- c(exons[hits.tx$subjectHits.x]$exon_id, exons[hits.tx$subjectHits.y]$exon_id)
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
        hits.start.idx <- stringr::str_detect(unstrsplit(exons[S4Vectors::subjectHits(hits.start)]$tx_name), 
                                              paste(tx.rank$tx_name, collapse = "|"))
        hits.end.idx <- stringr::str_detect(unstrsplit(exons[S4Vectors::subjectHits(hits.end)]$tx_name),
                                            paste(tx.rank$tx_name, collapse = "|"))
        
        
        
        insSite.gr <- c(gr[S4Vectors::queryHits(hits.start)[hits.start.idx]], 
                        partner(gr)[S4Vectors::queryHits(hits.end)[hits.end.idx]])
        insSite.gr$exon <- c(exons[S4Vectors::subjectHits(hits.start)[hits.start.idx]]$exon_id,
                             exons[S4Vectors::subjectHits(hits.end)[hits.end.idx]]$exon_id)
        insSite.gr$txs <- c(exons[S4Vectors::subjectHits(hits.start)[hits.start.idx]]$tx_name,
                            exons[S4Vectors::subjectHits(hits.end)[hits.end.idx]]$tx_name)
        #insSite.gr$txs <- insSite.gr$txs[sapply(insSite.gr$txs, function(x){x %in% tx.rank$tx_name})]line
        
        
        insSite.gr <- insSite.gr[!names(insSite.gr) %in% names(rt.gr)]
        insSite.gr <- c(insSite.gr, gr[names(partner(gr)) %in% names(insSite.gr)])
        
        #TODO: add L1/Alu annotation for insertion site filtering.
        
        
        return(GRangesList(insSite = insSite.gr, rt = rt.gr))
    }
}