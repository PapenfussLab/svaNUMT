# NA12878 analysis
library(ggplot)
library(stringr)

hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

gridss <- readVcf("W:/i.data.na12878/f3668a2e72832dbc3201eb771e8ec7a1.vcf", "hg19")
gridss <- gridss[abs(breakpointRanges(gridss)$svLen %na% 1000000000) > 50,]
ggr <- breakpointRanges(gridss)
ggr <- ggr[!str_detect(ggr$FILTER, "LOW_QUAL"),]

ghom <- referenceHomology(ggr, hg19)

ggplot(ghom, aes(x=exacthomlen, y=inexacthomlen)) +
	geom_point(aes(shape=ifelse(seqnames(ggr)==seqnames(partner(ggr)), "intra", "inter"))) +
	geom_density2d() +
	geom_jitter() +
	coord_cartesian(xlim=c(0,50), ylim=c(0,100)) +
	coord_cartesian(xlim=c(0,200), ylim=c(0,600)) +
	labs("Homology length")

ggplot(ghom, aes(x=inexacthomlen)) + geom_histogram(binwidth=1)
ggplot(ghom, aes(x=exacthomlen)) + geom_histogram(binwidth=1)

ggplot(ghom, aes(x=inexacthomlen, y=inexactscore)) +
	geom_point(aes(shape=ifelse(seqnames(ggr)==seqnames(partner(ggr)), "intra", "inter"))) +
	geom_density2d() +
	geom_jitter() +
	coord_cartesian(xlim=c(0,50), ylim=c(0,100)) +
	coord_cartesian(xlim=c(0,200), ylim=c(0,600)) +
	labs("Homology length")

nchom <- NULL
ncgr <- NULL
for(sample in c("W:/projects/liposarcoma/data/gridss/778/778.vcf",
	"W:/projects/liposarcoma/data/gridss/GOT3/bwa/GOT3-bwa.vcf",
	"W:/projects/liposarcoma/data/gridss/T1000/T1000")) {
	vcf <- readVcf(sample, "hg19")
	vcf <- vcf[abs(breakpointRanges(vcf)$svLen %na% 1000000000) > 50,]
	gr <- breakpointRanges(vcf)
	gr <- gr[!str_detect(gr$FILTER, "LOW_QUAL"),]
	hom <- referenceHomology(gr, hg19)
	gr$sample <- sample
	hom$sample <- sample
	ncgr <- c(ncgr, gr)
	nchom <- rbind(nchom, hom)
}
