.elementExtract.List <- function(x, offset=1) {
	lengths <- elementLengths(x)
	flat <- BiocGenerics::unlist(x)
	hasValue <- lengths >= offset
	flatOffset <- head(c(1, 1 + cumsum(lengths)), -1) + offset - 1
	flatOffset[!hasValue] <- length(flat) + 1 # out of bounds
	# need to strip XStringSet since that throws an error
	# on out of bounds instead of returning a correctly typed NA
	return(.unXStringSet(flat)[flatOffset])
}
.elementExtract.ANY <- function(x, offset=1) {
	if (is.null(x)) return(x)
	if (is.vector(x)) {
		if (offset==1) return(x)
		return(x[rep(length(x) + 1, length(x))])
	}
	result <- sapply(x, function(r) r[offset], USE.NAMES=FALSE)
	return(result)
}
.elementExtract.XStringSet <- function(x, offset=1) {
	return(.elementExtract.ANY(as.character(x), offset))
}
#' Extracts the element of each element at the given position
#'
#' @param x list-like object
#' @param offset offset of list
#' @export
setGeneric("elementExtract", function(x, offset=1) standardGeneric("elementExtract"))
setMethod("elementExtract", "XStringSet", .elementExtract.XStringSet)
setMethod("elementExtract", "List", .elementExtract.List)
setMethod("elementExtract", "ANY", .elementExtract.ANY)

#' converts an XStringSet to a character
setGeneric(".unXStringSet", function(x) x)
setMethod(".unXStringSet", "XStringSet", function(x) as.character(x))


#' Replaces the NA values in a with corresponding values in b
#' @export
'%na%' <- function(a, b) {
	if (is.null(a) || length(a) == 0) return(b)
	if (is.null(b) || length(b) == 0) return(a)
	return(ifelse(is.na(a), b, a))
}

#' Uses b if a is NULL
#' @export
'%null%' <- function(a, b) {
	if (is.null(a)) return(b)
	return (a)
}
