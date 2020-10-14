context('event detection functions')
na12878 <- readVcf(.testfile("gridss-na12878.vcf"), "")
numt <- readVcf(.testfile("MT.vcf"), "")
colo829 <- readVcf(system.file("extdata", "diploidSV.vcf", package = "StructuralVariantAnnotation"))

#NUMT detection
test_that("numtDetect returns nothing with no nu-mt fusion", {
    gr <- suppressWarnings(breakpointRanges(na12878))
    expect_equal(NULL,numtDetect(gr))
    expect_message(numtDetect(gr), 
                   "There is no NUMT event detected. Check whether 'chrM' or 'MT' is present in the VCF.")
})

test_that("numtDetect returns ", {
    gr <- numtDetect(breakpointRanges(numt))
    expect_equal(c("11", "11", "MT", "MT"), as.character(seqnames(gr)))
    expect_equal(c(49883572, 49883585, 16093, 16493), start(gr))
    expect_equal(c(49883572, 49883585, 16093, 16493), end(gr))
    expect_equal(c("+", "-", "-", "+"), as.character(strand(gr)))
    expect_equal(CharacterList("gridss9h", "gridss39h", "gridss9h", "gridss39h"), gr$candidatePartnerId)
})
#TODO: multiple candidatePartnerIds

#RT detection
#vcf <- VariantAnnotation::readVcf(vcf.file, "hg19")
#gr <- breakpointRanges(vcf, nominalPosition=TRUE)
#rt <- rtDetect(gr, genes, maxgap=30, minscore=0.6)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes <- TxDb.Hsapiens.UCSC.hg19.knownGene
test_that("RT detection returns two Granges", {
    gr <- breakpointRanges(colo829, nominalPosition=TRUE)
    rt <- rtDetect(gr, genes, maxgap=50, minscore=0.3)
    expect_equal(2,length(rt))
})

