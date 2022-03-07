context('event detection functions')
library(readr)
library(BSgenome.Hsapiens.UCSC.hg19)
na12878 <- readVcf(system.file("extdata", "gridss-na12878.vcf", package = "svaNUMT"))
numt <- readVcf(system.file("extdata", "MT.vcf", package = "svaNUMT"))
#known NUMTs as GRanges
numtS <- read_table(system.file("extdata", "numtS.txt", package = "svaNUMT"), col_names = FALSE)
colnames(numtS) <- c("bin", "seqnames", "start", "end", "name", "score", "strand")
numtS <- `seqlevelsStyle<-`(GRanges(numtS), "NCBI")
#mitochondria reference genome
genomeMT <- BSgenome.Hsapiens.UCSC.hg19$chrMT
#default numtDetect() parameters
max_ins_dist=10
maxgap_numtS=10 
min_len=20
min.Align=0.8

#NUMT detection
test_that("numtDetect returns an empty list and messages with no nu-mt fusion", {
    gr <- suppressWarnings(breakpointRanges(na12878))
    expect_equal(list(MT=NULL, known=NULL, insSeq=NULL),
                 numtDetect(gr, numtS, genomeMT, 
                                 max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                                 min_len=min_len, min.Align=min.Align))
    expect_message(numtDetect(gr, numtS, genomeMT, 
                              max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                              min_len=min_len, min.Align=min.Align), 
                   "There is no NUMT event detected. Check whether 'chrM' or 'MT' is present in the VCF.")
    expect_message(numtDetect(gr, numtS, genomeMT, 
                              max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                              min_len=min_len, min.Align=min.Align), 
                   "There is no MT sequence from known NUMT events detected.")
    expect_message(numtDetect(gr, numtS, genomeMT, 
                              max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                              min_len=min_len, min.Align=min.Align), 
                   "There is no NUMT event detected by insertion sequences.")
})

test_that("numtDetect returns an empty and message when no paired candidate insertion site is found within the given threshold", {
    gr <- suppressWarnings(breakpointRanges(numt))
    expect_equal(list(MT=NULL, known=NULL, insSeq=NULL), 
                 numtDetect(gr, numtS, genomeMT, 
                                  max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                                  min_len=min_len, min.Align=min.Align))
    expect_message(numtDetect(gr, numtS, genomeMT, 
                              max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                              min_len=min_len, min.Align=min.Align), 
                   "There is no NUMT event detected. Paired candidate insertion site not found.")
    expect_message(numtDetect(gr, numtS, genomeMT, 
                              max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                              min_len=min_len, min.Align=min.Align), 
                   "There is no MT sequence from known NUMT events detected.")
    expect_message(numtDetect(gr, numtS, genomeMT, 
                              max_ins_dist = 10, maxgap_numtS=maxgap_numtS, 
                              min_len=min_len, min.Align=min.Align), 
                   "There is no NUMT event detected by insertion sequences.")
    # expect_equal(c("11", "11", "MT", "MT"), as.character(seqnames(gr)))
    # expect_equal(c(49883572, 49883585, 16093, 16493), start(gr))
    # expect_equal(c(49883572, 49883585, 16093, 16493), end(gr))
    # expect_equal(c("+", "-", "-", "+"), as.character(strand(gr)))
    # expect_equal(CharacterList("gridss9h", "gridss39h", "gridss9h", "gridss39h"), gr$candidatePartnerId)
})

test_that("numtDetect returns a list element MT which consists of 2 lists when 
          candidate NUMTs by MT breakpoints are found", {
    gr <- suppressWarnings(breakpointRanges(numt))
    NUMT <- numtDetect(gr, numtS, genomeMT, 
               max_ins_dist = 20, maxgap_numtS=maxgap_numtS, 
               min_len=min_len, min.Align=min.Align)
    expect_equal(2, length(NUMT$MT))
})

