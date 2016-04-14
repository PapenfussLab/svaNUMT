#' Expands any columns with multiple values (such as Lists) into separate columns
#'
#' @export
setGeneric("unpack",
		   function(x, ...)
		   	standardGeneric("unpack")
)
setMethod("unpack", "DataFrame",
		  function(x, ...)
		  	.unpack.DataFrame(x)
)
setMethod("unpack", "data.frame",
  	  function(x, ...)
  	  	.unpack.DataFrame(x)
)
setMethod("as.matrix", "List",
		  function(x, ...)
		  	.as.matrix.List(x)
)

.as.matrix.List <- function(x) {
	ncolx <- max(lengths(x))
	nrowx <- length(x)
	occupiedcells <- inverse.rle(list(lengths=c(rbind(lengths(x), ncolx - lengths(x))), values=rep(c(TRUE, FALSE), length.out=2*nrowx)))
	data <- rep(NA, ncolx * nrowx)
	data[occupiedcells] <- unlist(x)
	mat <- matrix(data, nrow=nrowx, ncol=ncolx, byrow=TRUE)
	return(mat)
}

#' Returns a list of aggegration functions for the given
#' variant caller
aggregatesFor <- function(caller) {

}
#' Unpacks the columns of the given data frame
#'
unpack.data.frame(x, caller, F=aggregatesFor(caller)) {

}
.unpack.DataFrame <- function(x, ...) {
	unpack(as.data.frame(x), ...)
}
.unpack.data.frame <- function(x, f=list(columnName=F)) {un
	for (cname in names(x)) {
		cx <- x[cname]
		if (is.list(cx)) {
			mat <-
		}
	}
}

