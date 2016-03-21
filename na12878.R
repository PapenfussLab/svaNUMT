# NA12878 analysis
library(ggplot2)
library(stringr)
library(VariantAnnotation)

hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

gridss <- readVcf("~/i/data.na12878/f3668a2e72832dbc3201eb771e8ec7a1.vcf", "hg19")
gridss <- gridss[abs(breakpointRanges(gridss)$svLen %na% 1000000000) > 50,]
ggr <- breakpointRanges(gridss)
ggr <- ggr[!str_detect(ggr$FILTER, "LOW_QUAL") & ggr$QUAL >= 300,]
ghom <- referenceHomology(ggr, hg19)
ghom$sample <- "NA12878"
ggr$sample <- "NA12878"


nchom <- NULL
ncgr <- NULL
#setwd("W:/projects/liposarcoma/data/gridss")
for(sample in c("W:/projects/liposarcoma/data/gridss/778/778.vcf",
	"W:/projects/liposarcoma/data/gridss/GOT3/bwa/GOT3-bwa.vcf",
	"W:/projects/liposarcoma/data/gridss/T1000/T1000.vcf")) {
	vcf <- readVcf(sample, "hg19")
	vcf <- vcf[abs(breakpointRanges(vcf)$svLen %na% 1000000000) > 50,]
	gr <- breakpointRanges(vcf)
	gr <- gr[!str_detect(gr$FILTER, "LOW_QUAL"),]
	gr <- gr[gr$QUAL >= 300,]
	hom <- referenceHomology(gr, hg19)
	gr$sample <- sample
	hom$sample <- sample
	ncgr <- c(ncgr, gr)
	nchom <- rbind(nchom, hom)
}


gr <- c(ncgr, ggr)
hom <- rbind(nchom, ghom)
hom$loc <- ifelse(seqnames(gr)==seqnames(partner(gr)), "intra", "inter")
ggplot(hom, aes(x=exacthomlen, y=inexacthomlen)) +
	geom_point(aes(color=sample, shape=loc) +
	geom_density2d() +
	#geom_jitter() +
	coord_cartesian(xlim=c(0,50), ylim=c(0,100)) +
	coord_cartesian(xlim=c(0,200), ylim=c(0,600)) +
	labs("Homology length")

ggplot(hom, aes(x=inexacthomlen, color=sample)) +
	facet_wrap(~sample) +
	geom_histogram(binwidth=1) +
	coord_cartesian(xlim=c(0,25))
ggplot(hom, aes(x=exacthomlen, color=sample)) +
	geom_density() +
	geom_histogram(binwidth=1) +
	coord_cartesian(xlim=c(0,100))

ggplot(hom, aes(x=inexactscore / inexacthomlen, y=inexacthomlen, color=loc)) +
	geom_point() +
	geom_rug(col=rgb(.5,0,0,alpha=.05), sides="l") +
	facet_wrap(~sample) +
	coord_cartesian(xlim=c(0,2), ylim=c(0,600)) +
	labs(title="Inexact homology length", y="inexact homlen", x="alignment score per base (match=2, mismatch=-6)")


