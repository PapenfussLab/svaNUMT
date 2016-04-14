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
setMethod("unpack", "VCF",
		  function(x, ...)
		  	.unpack.DataFrame(info(x))
)

.as.matrix.list <- function(x) {
	ncolx <- max(lengths(x))
	nrowx <- length(x)
	occupiedcells <- inverse.rle(list(lengths=c(rbind(lengths(x), ncolx - lengths(x))), values=rep(c(TRUE, FALSE), length.out=2*nrowx)))
	data <- rep(NA, ncolx * nrowx)
	data[occupiedcells] <- unlist(x)
	mat <- matrix(data, nrow=nrowx, ncol=ncolx, byrow=TRUE)
	return(mat)
}
.aggregate_functions = list(
	"gridss" = list(
		ASCRP=sum,
		ASCRR=sum,
		ASRP=sum,
		ASSR=sum,
		BSC=sum,
		BSCQ=sum,
		BUM=sum,
		BUMQ=sum,
		#HOMLEN=max,
		#HOMSEQ=sum,
		REF=sum,
		REFPAIR=sum,
		RP=sum,
		RPQ=sum,
		RSR=sum,
		RSRQ=sum,
		SR=sum,
		SRQ=sum
	)
)
#' Returns a list of aggegration functions for the given
#' variant caller
#' @export
aggregatesFor <- function(caller) {
	.aggregate_functions[[caller]]
}
#' Unpacks the columns of the given data frame
#'
#'@export
calculate_aggregates <- function(x, caller="gridss", F=aggregatesFor(caller)) {

}
.unpack.DataFrame <- function(x, caller="gridss", F=aggregatesFor(caller), ...) {
	xdf <- as.data.frame(x)
	for (cname in names(x)) {
		cx <- x[cname]
		if (is.list(cx) || is(cx, "List")) {
			mat <- .as.matrix.list(cx)
			cdf <- as.data.frame(mat)
			names(cdf) <- rep(cname, length(names(cdf)))
			if (!is.null(F[[cname]]) {
				xdf[[cname]] <- F[[cname]](mat)
			}
			xdf <- cbind(xdf, cdf)
		}
	}
}

