context('event detection functions')
na12878 <- readVcf(system.file("extdata", "gridss-na12878.vcf", package = "NUMTDetect"))
numt <- readVcf(system.file("extdata", "MT.vcf", package = "NUMTDetect"))
#colo829 <- readVcf(system.file("extdata", "diploidSV.vcf", package = "NUMTDetect"))

#NUMT detection
test_that("numtDetect returns nothing with no nu-mt fusion", {
    gr <- suppressWarnings(breakpointRanges(na12878))
    expect_equal(NULL,numtDetect(gr))
    expect_message(numtDetect(gr), 
                   "There is no NUMT event detected. Check whether 'chrM' or 'MT' is present in the VCF.")
})

test_that("numtDetect returns nothing but message when no paired candidate insertion site is found within the given threshold", {
    gr <- suppressWarnings(breakpointRanges(numt))
    expect_equal(NULL, numtDetect(gr, 10))
    expect_message(numtDetect(gr, 10), 
                   "There is no NUMT event detected. Paired candidate insertion site not found.")
    
    # expect_equal(c("11", "11", "MT", "MT"), as.character(seqnames(gr)))
    # expect_equal(c(49883572, 49883585, 16093, 16493), start(gr))
    # expect_equal(c(49883572, 49883585, 16093, 16493), end(gr))
    # expect_equal(c("+", "-", "-", "+"), as.character(strand(gr)))
    # expect_equal(CharacterList("gridss9h", "gridss39h", "gridss9h", "gridss39h"), gr$candidatePartnerId)
})

test_that("numtDetect returns a list of 2 when candidate NUMTs are found", {
    gr <- suppressWarnings(breakpointRanges(numt))
    expect_equal(2, length(numtDetect(gr, 20)))
})

