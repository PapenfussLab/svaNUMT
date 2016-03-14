#' Extracts the element at the given position
.elementExtract <- function(x, offset=1) {
    #assert_that(is(x, "List"))
    result <- sapply(x, function(r) r[offset], USE.NAMES=FALSE)
    return(result)
}