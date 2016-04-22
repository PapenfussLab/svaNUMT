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
	rlens <- BiocGenerics::lengths(x)
	ncolx <- max(rlens)
	nrowx <- length(x)
	occupiedcells <- inverse.rle(list(lengths=c(rbind(rlens, ncolx - rlens)), values=rep(c(TRUE, FALSE), length.out=2*nrowx)))
	occupiedData <- BiocGenerics::unlist(x)
	# populate with NAs of the correct type
	data <- append(occupiedData[c()], rep(NA, ncolx * nrowx))
	data[occupiedcells] <- occupiedData
	mat <- matrix(data, nrow=nrowx, ncol=ncolx, byrow=TRUE)
	return(mat)
}
.drop <- function(x) { NULL }
.first <- function(x) { elementExtract(x) }
.zeroFill <- function(x) { x %na% 0 }
.emptyFill <- function(x) { x %na% "" }
.predefined_aggregates = list(
	"all" = list(
		SVLEN=.first,
		HOMLEN=.first,
		HOMSEQ=.first,
		CIPOS=.drop,
		CIEND=.drop
	),
	"gridss" = list(
		ASCRP=sum,
		ASCSR=sum,
		ASRP=sum,
		ASSR=sum,
		BSC=sum,
		BSCQ=sum,
		BUM=sum,
		BUMQ=sum,
		REF=sum,
		REFPAIR=sum,
		RP=sum,
		RPQ=sum,
		RSR=sum,
		RSRQ=sum,
		SR=sum,
		SRQ=sum,
		TEST=.drop,
		IHOMPOS=.drop,
		CIRPOS=.drop,
		RSI=.drop,
		SI=.drop,
		MATEID=.first
	)
)
.predefined_transforms = list(
	"all" = list(
		SVLEN=.first,
		HOMLEN=.zeroFill,
		HOMSEQ=.first,
		CIPOS=.zeroFill,
		CIEND=.zeroFill
	),
	"gridss" = list(
		# per sample
		ASCRP=.zeroFill,
		ASCSR=.zeroFill,
		ASRP=.zeroFill,
		ASSR=.zeroFill,
		BSC=.zeroFill,
		BSCQ=.zeroFill,
		BUM=.zeroFill,
		BUMQ=.zeroFill,
		HOMLEN=.zeroFill,
		HOMSEQ=function(x) { elementExtract(x) %na% "" },
		REF=.zeroFill,
		REFPAIR=.zeroFill,
		RP=.zeroFill,
		RPQ=.zeroFill,
		RSR=.zeroFill,
		RSRQ=.zeroFill,
		SR=.zeroFill,
		SRQ=.zeroFill,
		# single value fields
		AS=.zeroFill,
		ASQ=.zeroFill,
		BA=.zeroFill,
		BAQ=.zeroFill,
		BQ=.zeroFill,
		CQ=.zeroFill,
		RAS=.zeroFill,
		RASQ=.zeroFill,
		SC=.zeroFill,
		# interval fields
		IHOMPOS=.zeroFill,
		CIRPOS=.zeroFill
	)
)
#' Returns a list of aggegration functions for the given variant caller
#' @export
aggregateFunctionsFor <- function(caller) {
	c(.predefined_aggregates[["all"]], .predefined_aggregates[[caller]])
}
#' Returns a list of transformation functions for the given variant caller
#' @export
transformFunctionsFor <- function(caller) {
	c(.predefined_transforms[["all"]], .predefined_transforms[[caller]])
}
.unpack.DataFrame <- function(x, caller="gridss", transformFUN=transformFunctionsFor(caller), aggregateFUN=aggregateFunctionsFor(caller), ...) {
	xdf <- as.data.frame(x)
	for (cname in names(x)) {
		cx <- x[[cname]]
		if (is.list(cx) || is(cx, "List")) {
			# apply aggregation to list columns
			if (!is.null(aggregateFUN[[cname]])) {
				xdf[[cname]] <- aggregateFUN[[cname]](cx)
			}
			# split out list columns
			mat <- .as.matrix.list(cx)
			cdf <- as.data.frame(mat)
			if (!is.null(transformFUN[[cname]])) {
				for (i in seq_along(cdf)) {
					# apply transform to each unpacked column
					cdf[[i]] <- transformFUN[[cname]](cdf[[i]])
				}
			}
			names(cdf) <- rep(cname, length(cdf))
			xdf <- cbind(xdf, cdf)
		}
	}
	names(xdf) <- make.unique(names(xdf))
	for (cname in names(x)) {
		# apply transform to named column
		if (!is.null(transformFUN[[cname]]) && !is.null(xdf[[cname]])) {
			xdf[[cname]] <- transformFUN[[cname]](xdf[[cname]])
		}
	}
	return(xdf)
}

