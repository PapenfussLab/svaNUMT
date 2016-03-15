#' Extracts the element at the given position
#' TODO: turn this into a generic
#' Biostrings::XStringSetList
#' S4Vectors::List
.elementExtract <- function(x, offset=1) {
	if (is(x, "DNAStringSetList")) {
	} else if ()
    result <- sapply(x, function(r) r[offset], USE.NAMES=FALSE)
    return(result)
}
