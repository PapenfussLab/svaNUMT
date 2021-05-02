#' Calculating MT sequence length.
#'
#' @details
#' This function calculate the length of MT sequence length with BND notations.
#' @param bnd.start starting breakend of the MT sequence.
#' @param bnd.end ending breakend of the MT sequence.
#' @param chrM.len length of the reference MT genome.
#' @keywords internal
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