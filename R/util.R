.elementExtract.List <- function(x, offset=1) {
	lengths <- S4Vectors::elementNROWS(x)
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
setGeneric("elementExtract", function(x, offset=1) standardGeneric("elementExtract"))
setMethod("elementExtract", "XStringSet", .elementExtract.XStringSet)
setMethod("elementExtract", "List", .elementExtract.List)
setMethod("elementExtract", "ANY", .elementExtract.ANY)

#' converts an XStringSet to a character
#' @param x an XStringSet.
setGeneric(".unXStringSet", function(x) x)
setMethod(".unXStringSet", "XStringSet", function(x) as.character(x))


#' Replaces the NA values in a with corresponding values in b
#' @param a,b objects to be tested or coerced.
'%na%' <- function(a, b) {
	if (is.null(a) || length(a) == 0) return(b)
	if (is.null(b) || length(b) == 0) return(a)
	return(ifelse(is.na(a), b, a))
}

#' Uses b if a is NULL
#'  @param a,b objects to be tested or coerced.
'%null%' <- function(a, b) {
	if (is.null(a)) return(b)
	return (a)
}

#' vectorised pairwise longest common prefix
#' Returns the length of the longest common prefix for
#' each string pair
#' @param s1,s2 A pair of strings.
#' @param ignore.case Whether cases in the strings should be ignored.
.pairwiseLCPrefix <- function(s1, s2, ignore.case=FALSE) {
	s1 <- as.character(s1)
	s2 <- as.character(s2)
	if (ignore.case) {
		s1 <- toupper(s1)
		s2 <- toupper(s2)
	}
	prefixLength <- rep(0, max(length(s1), length(s2)))
	matchi <- TRUE
	i <- 1
	while (any(matchi)) {
		s1i <- substring(s1, i, i)
		s2i <- substring(s2, i, i)
		matchi <- s1i != "" & s1i == s2i
		prefixLength <- prefixLength + as.integer(matchi)
		i <- i + 1
	}
	return(prefixLength)
}




