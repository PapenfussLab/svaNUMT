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
svqsc_train <- function(traininggr, trainingdf, truthgr, considerDuplicateCallsTrue=FALSE, requiredSupportingReads=1, allowsPartialHits=FALSE, ...) {
	assert_that(!any(duplicated(names(trainingdf))))
	truth <- .svqsc_matches_truth(traininggr, trainingdf, truthgr, considerDuplicateCallsTrue, requiredSupportingReads, allowsPartialHits, ...)
	model <- .svqsc_generate_model(trainingdf, truth)
	return(model)
}
.svqsc_matches_truth <- function(traininggr, trainingdf, truthgr, considerDuplicateCallsTrue, requiredSupportingReads, allowsPartialHits, ...) {
	hitscounts <- .svqsc_long_read_hits(traininggr, truthgr, considerDuplicateCallsTrue, ...)
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
	result <- rep(FALSE, times=nrow(trainingdf))
	names(result) <- row.names(trainingdf)
	if (nrow(hitdf) > 0) {
		result[traininggr[hitdf$name]$vcfId] <- TRUE
	}
	names(result) <- NULL
	assert_that(length(result) == nrow(trainingdf))
	return(result)
}
#' Calls variants
#' @param considerDuplicateCallsTrue count each read towards all matching calls.
.svqsc_long_read_hits <- function(callgr, truthgr, considerDuplicateCallsTrue, ...) {
	hitscounts <- rep(0, length(callgr))
	hits <- findBreakpointOverlaps(callgr, truthgr, ...)
	if (considerDuplicateCallsTrue) {
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
	#assert_that(!any(is.na(glmmat)))
	# Missing values default to zero
	glmmat[is.na(glmmat)] <- 0
	#matrix(as.numeric(unlist(df)),nrow=nrow(df))
	#colnames(glmmat) <- names(df)
	return(glmmat)
}
#' generates a calibrated model from a numeric data frame with a tp field
.svqsc_generate_model <- function(trainingdf, truth) {
	assert_that(nrow(trainingdf) == length(truth))
	# subset so tp and fp sets are of similar size
	trainingdf$tp <- truth
	subsetdf <- rbind(
		trainingdf %>% filter(tp),
		trainingdf %>% filter(!tp) %>% sample_n(max(100,sum(trainingdf$tp))))

	glmmat <- .svqsc_df2matrix(subsetdf %>% select(-tp))
	glmmod <- glmnet::cv.glmnet(glmmat, y=subsetdf$tp, alpha=1, family='binomial')
	# TODO: is this calibrated?
	return (glmmod)
}
#' Score the given variants according to the given model
#' @param df breakpoint annotations
#' @param model scoring model generated from \code{\link{svqsc_train}}
svqsc_score <- function(model, df) {
	pred <- 1 - predict(model, newx=as.matrix(df), type="response", s="lambda.1se")
	phredpred <- -10*log10(pred)
	return(phredpred)
}

#' Generates a precision recall plot for QUAL, and for the model
.svqsc_precision_recall <- function(vcf, truthgr, dftransform, requiredSupportingReads, considerDuplicateCallsTrue=FALSE, allowsPartialHits=FALSE, intrachromosomalOnly=TRUE) {
	# filters: 50bp
	vcf <- vcf[is.na(StructuralVariantAnnotation:::.svLen(vcf)) | abs(StructuralVariantAnnotation:::.svLen(vcf)) >= minsize,]
	gr <- breakpointRanges(vcf)
	if (intrachromosomalOnly) {
		intravcfids <- gr$vcfId[as.logical(seqnames(gr) == seqnames(partner(gr)))]
		intravcfids <- intravcfids[!duplicated(intravcfids)]
		vcf <- vcf[intravcfids,]
		gr <- breakpointRanges(vcf)
	}
	if (any(is.na(gr$QUAL))) {
		stop("Missing QUAL score for at least one breakpoint")
	}
	df <- StructuralVariantAnnotation::unpack(vcf)
	df$QUAL <- 0
	df[gr$vcfId,]$QUAL <- gr$QUAL
	dfQUAL <- df$QUAL
	df <- dftransform(df, gr)
	truth <- .svqsc_matches_truth(gr, df, truthgr, considerDuplicateCallsTrue, requiredSupportingReads, allowsPartialHits, maxgap=maxgap, ignore.strand=ignore.strand)
	model <- svqsc_train(gr, df, truthgr, maxgap=maxgap, ignore.strand=ignore.strand, requiredSupportingReads)
	modelpred <- svqsc_score(model, df)
	#pred <- df %>%
	#	dplyr::mutate(estimator="QUAL") %>%
	#	dplyr::select(estimator, QUAL, tp) %>%
	#	rbind(data.frame(
	#		estimator=rep(colnames(allpred), each=nrow(allpred)),
	#		QUAL=as.vector(allpred),
	#		tp=rep(df$tp, times=ncol(allpred))))
	qual_score <- data.frame(estimator="QUAL", QUAL=dfQUAL, tp=truth)
	svqsc_score <- data.frame(estimator="lambda.1se", QUAL=unname(modelpred), tp=truth)
	pred <- rbind(qual_score, svqsc_score)
	roc <- pred %>%
		dplyr::group_by(estimator) %>%
		dplyr::arrange(desc(QUAL)) %>%
		dplyr::mutate(fp=!tp) %>%
		dplyr::mutate(tp=cumsum(tp), fp=cumsum(fp)) %>%
		dplyr::group_by(estimator, QUAL) %>%
		dplyr::summarise(tp=max(tp), fp=max(fp)) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(precision=tp / (tp + fp), n=tp + fp)

	roc <- roc %>%
	  dplyr::group_by(estimator) %>%
	  dplyr::arrange(n) %>%
	  dplyr::filter(
	    # keep start/end
	    is.na(dplyr::lag(tp)) | is.na(dplyr::lead(tp)) |
	      # keep group transitions (TODO: is there a way to make lead/lag across group_by return NA?)
	      estimator != dplyr::lag(estimator) |
	      estimator != dplyr::lead(estimator) |
	      # slopes not equal dx1/dy1 != dx2/dy2 -> dx1*dy2 != dx2*dy1
	      (tp - dplyr::lag(tp))*(dplyr::lead(fp) - dplyr::lag(fp)) != (dplyr::lead(tp) - dplyr::lag(tp))*(fp - dplyr::lag(fp)) |
	      # less than 10 calls wide
	      dplyr::lead(tp) - dplyr::lag(tp) > 10 |
	      # keep every 5th row
	      row_number() %% 10 == 0)
	# lossy removal of points with least change (shinyCache.R)
	for (k in c(4, 16, 32, 64)) {
		roc <- roc %>%
			dplyr::group_by(estimator) %>%
			dplyr::arrange(n) %>%
			dplyr::filter(
				is.na(dplyr::lag(tp)) | is.na(dplyr::lead(tp)) |
				estimator != dplyr::lag(estimator) |
				estimator != dplyr::lead(estimator) |
				# remove points with least amount of change
				dplyr::lead(tp) - dplyr::lag(tp) + dplyr::lead(fp) - dplyr::lag(fp) > k |
				# keep every 5th to prevent removal of large segments
				row_number() %% 5 == 0
			) %>%
			dplyr::ungroup()
	}
	precisionRecallPlot <- ggplot(roc) +
		aes(y=precision, x=tp, color=estimator) +
		geom_line() +
		geom_line(data=roc %>% filter(estimator=="QUAL"), size=2)
	return(list(model=model, roc=roc, precisionRecallPlot=precisionRecallPlot))
}


