#' Generates a calibrated variant scoring model
#'
#' @param traininggr breakpoints used to train the model
#' @param trainingdf annotation used for training
#' @param truthgr breakpoints considered true positives.
#' @param requiredSupportingReads truthgr support count required to be called as a true positive.
#' When comparing against a known truth call set, use 1.
#' When comparing directly against long read split alignments, a value greater than 1 is appropriate.
#' @param allowsPartialHits count truth matches at the event level.
#' If TRUE, a variant is called as a true positive if supported by requiredSupportingReads in truthgr
#' If FALSE, a variant is called as a true positive if all constituent breakpoints are supported by requiredSupportingReads in truthgr
#' @return Variant scoring model for use in \code{\link{svqsc_score}}
#' @export
svqsc_train <- function(traininggr, trainingdf, truthgr, countAllMatchingCalls=FALSE, requiredSupportingReads=1, allowsPartialHits=FALSE, ...) {
	dftp <- .svqsc_annotate_tp(traininggr, trainingdf, truthgr, countAllMatchingCalls, requiredSupportingReads, allowsPartialHits, ...)
	return(svqsc_generate_model(dftp))
}
.svqsc_annotate_tp <- function(traininggr, trainingdf, truthgr, countAllMatchingCalls, requiredSupportingReads, allowsPartialHits, ...) {
	hitscounts <- .svqsc_long_read_hits(traininggr, truthgr, countAllMatchingCalls, ...)
	if (allowsPartialHits) {
		hitdf <- data.frame(
			name=names(traininggr),
			hitscounts=hitscounts) %>%
			dplyr::group_by(name) %>%
			dplyr::summarise(hitscounts=sum(hitscounts)) %>%
			dplyr::mutate(tp=hitscounts >= requiredSupportingReads) %>%
			dplyr::filter(tp)
	} else {
		hitdf <- data.frame(
			name=names(traininggr),
			hitscounts=hitscounts,
			tp=hitscounts >= requiredSupportingReads) %>%
			dplyr::group_by(name) %>%
			dplyr::summarise(tp=all(tp), hitscounts=sum(hitscounts)) %>%
			dplyr::filter(tp)
	}
	trainingdf$tp <- FALSE
	trainingdf[traininggr[hitdf$name]$vcfId,]$tp <- TRUE
	return(trainingdf)
}
#' Calls variants
#' @param countAllMatchingCalls count each read towards all matching calls.
.svqsc_long_read_hits <- function(callgr, truthgr, countAllMatchingCalls=FALSE, ...) {
	hitscounts <- rep(0, length(callgr))
	hits <- findBreakpointOverlaps(callgr, truthgr, ...)
	if (countAllMatchingCalls) {
		hits <- hits %>%
	      dplyr::group_by(queryHits) %>%
	      dplyr::summarise(n=n())
	} else {
		# assign supporting evidence to the call with the highest QUAL
		hits$QUAL <- callgr$QUAL[hits$queryHits]
	    hits <- hits %>%
	      dplyr::arrange(desc(QUAL), queryHits) %>%
	      dplyr::distinct(subjectHits, .keep_all=TRUE) %>%
	      dplyr::group_by(queryHits) %>%
	      dplyr::summarise(n=n())
	}
    hitscounts[hits$queryHits] <- hits$n
    return(hitscounts)
}
#' Convert the annotation data frame to a matrix
#' TODO: how do we generically handle non-numeric values? Do we just drop them?
.svqsc_df2matrix <- function(df) {
	glmmat <- as.matrix(df)
	#matrix(as.numeric(unlist(df)),nrow=nrow(df))
	#colnames(glmmat) <- names(df)
	return(glmmat)
}
#' generates a calibrated model from a numeric data frame with a tp field
.svqsc_generate_model <- function(trainingdf) {
	glmmat <- .svqsc_df2matrix(trainingdf %>% dplyr::select(-tp))
	glmmod <- glmnet::glmnet(glmmat, y=trainingdf$tp, alpha=1, family='binomial')
	# TODO: is this calibrated?
	return (glmmod)
}
#' Score the given variants according to the given model
#' @param df breakpoint annotations
#' @param model scoring model generated from \code{\link{svqsc_train}}
svqsc_score <- function(model, df, ...) {
	allpred <- predict(glmmod, newx=as.matrix(df), type="response", ...)
	return(allpred)
}





