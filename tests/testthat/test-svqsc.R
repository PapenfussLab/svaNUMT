library(dplyr)
library(ggplot2)

minsize <- 51
maxgap <- 200
ignore.strand <- TRUE
sizemargin <- 0.25
requiredSupportingReads <- 5
#truthgr <- bedpe2breakpointgr("C:/dev/sv_benchmark/input.na12878/lumpy-Mills2012-call-set.bedpe")
truthgr <- bedpe2breakpointgr("C:/dev/input.na12878/longread/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.noid.sv.bedpe.gz")
wgEncodeDacMapabilityConsensusExcludable <- import("C:/dev/sv_benchmark/input.na12878/wgEncodeDacMapabilityConsensusExcludable.bed")
seqlevelsStyle(wgEncodeDacMapabilityConsensusExcludable) <- "UCSC"
seqlevelsStyle(truthgr) <- "UCSC"

encodeAnnotate <- function(df, gr) {
	df$wgEncodeDacMapabilityConsensusExcludable <- FALSE
	df[gr$vcfId[overlapsAny(gr, truthgr)],]$wgEncodeDacMapabilityConsensusExcludable <- TRUE
	return(df)
}
go <- function(vcf, dftransform) {
	return(.svqsc_precision_recall(vcf, truthgr, requiredSupportingReads, dftransform=function(df, gr) return(encodeAnnotate(dftransform(df, gr), gr))))
}
result <- list()

breakdancervcf <- readVcf("C:/dev/svqsc/breakdancer.vcf.gz", "hg19")
info(breakdancervcf)[!is.na(info(breakdancervcf)$EVENT)]$Chr1

result[["breakdancer"]] <- go(breakdancervcf, dftransform=function(df, gr) df %>%
		dplyr::select(IMPRECISE, SVLEN, Size, Score, num_Reads))

dellyvcf <- readVcf("C:/dev/svqsc/delly.vcf.gz", "hg19")
fixed(dellyvcf)$QUAL <- (info(dellyvcf)$PE %na% 0) + (info(dellyvcf)$SR %na% 0)
result[["delly"]] <- go(dellyvcf, dftransform=function(df, gr) df %>%
		dplyr::select(PE, MAPQ, SR, SRQ, IMPRECISE, CIPOS))

socratesvcf <- readVcf("C:/dev/svqsc/socrates.vcf.gz", "hg19")
result[["socrates"]] <- go(socratesvcf, dftransform=function(df, gr) df %>%
		# todo mutate(ANCHCONSLEN=nchar(ANCHCONS), REALNCONSLEN=nchar(REALNCONS))
		dplyr::select(QUAL, NLSC, NSSC, BLSC, BSSC, LSSC))






# debug
remove(vcf, df, gr)
considerDuplicateCallsTrue <- FALSE
requiredSupportingReads <- 3
allowsPartialHits <- TRUE
intrachromosomalOnly <- TRUE
traininggr <- gr
trainingdf <- df
hitscounts <- .svqsc_long_read_hits(traininggr, truthgr, considerDuplicateCallsTrue, maxgap=maxgap, ignore.strand=ignore.strand)
trainingdf <- .svqsc_annotate_tp(traininggr, trainingdf, truthgr, considerDuplicateCallsTrue, requiredSupportingReads, allowsPartialHits, maxgap=maxgap, ignore.strand=ignore.strand)
vcf <- socratesvcf


cv <- cv.glmnet(as.matrix(trainingdf %>% select(-tp)), trainingdf$tp, alpha=1, family='binomial')
pred <- predict(cv, newx=as.matrix(df), type="response", s="lambda.1se")



library(caret) #install.packages('caret', dependencies=TRUE)
x <- as.matrix(trainingdf %>% select(-tp))
y <- trainingdf$tp
trc = trainControl(method="cv", number=10)
fitM = train(x, y, trControl=trc, method="glmnet", family="binomial", metric="Accuracy")



