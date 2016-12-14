maxgap <- 200
ignore.strand <- TRUE
sizemargin <- 0.25
minlongreads <- 1
truthgr <- bedpe2breakpointgr("C:/dev/sv_benchmark/input.na12878/lumpy-Mills2012-call-set.bedpe")
dellyvcf <- readVcf("C:/dev/data.na12878/f37c9c1cc36a2297146c8e309bdab35d.vcf", "hg19")



vcf <- dellyvcf
# filters: 50bp
vcf <- vcf[abs(StructuralVariantAnnotation:::.svLen(vcf)) > 50,]
# TODO: remove interchromosomal
gr <- breakpointRanges(dellyvcf)
df <- StructuralVariantAnnotation::unpack(dellyvcf) %>%
	dplyr::select(PE, MAPQ, SR, SRQ, IMPRECISE, CIPOS)

model <- svqsc_train(gr, df, truthgr, maxgap=maxgap, ignore.strand=ignore.strand, requiredSupportingReads=minlongreads)
allpred <- svqsc_score(model, df)
library(dplyr)
library(ggplot2)
df <- .svqsc_annotate_tp(gr, df, truthgr, countAllMatchingCalls, requiredSupportingReads, allowsPartialHits, maxgap=maxgap, ignore.strand=ignore.strand)
df$QUAL <- df$PE + df$SR

pred <- df %>%
	mutate(estimator="QUAL") %>%
	select(estimator, QUAL, tp) %>%
	rbind(data.frame(
		estimator=rep(colnames(allpred), each=nrow(allpred)),
		QUAL=as.vector(allpred),
		tp=rep(df$tp, times=ncol(allpred))))
roc <- pred %>%
	group_by(estimator) %>%
	arrange(desc(QUAL)) %>%
	mutate(fp=!tp) %>%
	mutate(tp=cumsum(tp), fp=cumsum(fp)) %>%
	group_by(estimator, QUAL) %>%
	summarise(tp=max(tp), fp=max(fp)) %>%
	ungroup() %>%
	mutate(precision=tp / (tp + fp), n=tp + fp)

roc <- roc %>%
  group_by(estimator) %>%
  arrange(n) %>%
  filter(
    # keep start/end
    is.na(lag(tp)) | is.na(lead(tp)) |
      # keep group transitions (TODO: is there a way to make lead/lag across group_by return NA?)
      estimator != lag(estimator) |
      estimator != lead(estimator) |
      # slopes not equal dx1/dy1 != dx2/dy2 -> dx1*dy2 != dx2*dy1
      (tp - lag(tp))*(lead(fp) - lag(fp)) != (lead(tp) - lag(tp))*(fp - lag(fp)) |
      # less than 10 calls wide
      lead(tp) - lag(tp) > 10 |
      # keep every 5th row
      row_number() %% 10 == 0)
# lossy removal of points with least change (shinyCache.R)
for (k in c(4, 16, 32, 64)) {
	roc <- roc %>%
		group_by(estimator) %>%
		arrange(n) %>%
		filter(
			is.na(lag(tp)) | is.na(lead(tp)) |
			estimator != lag(estimator) |
			estimator != lead(estimator) |
			# remove points with least amount of change
			lead(tp) - lag(tp) + lead(fp) - lag(fp) > k |
			# keep every 5th to prevent removal of large segments
			row_number() %% 5 == 0
		) %>%
		ungroup()
}

ggplot(roc) +
	aes(y=precision, x=tp, color=estimator) +
	geom_line() +
	geom_line(data=roc %>% filter(estimator=="QUAL"), size=2)




# debug
countAllMatchingCalls <- FALSE
requiredSupportingReads <- 1
allowsPartialHits <- TRUE
traininggr <- gr
trainingdf <- df
callgr <- gr
hitscounts <- .svqsc_long_read_hits(traininggr, truthgr, countAllMatchingCalls, maxgap=maxgap, ignore.strand=ignore.strand)
trainingdf <- .svqsc_annotate_tp(traininggr, trainingdf, truthgr, countAllMatchingCalls, requiredSupportingReads, allowsPartialHits, maxgap=maxgap, ignore.strand=ignore.strand)

