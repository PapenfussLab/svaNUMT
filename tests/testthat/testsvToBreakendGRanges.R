
example <- readVcf(.testfile("vcf4.2.example.sv.vcf"), "")
simple <- readVcf(.testfile("simple.vcf"), "")

breakdancer <- readVcf(.testfile("breakdancer-1.4.5.vcf"), "")
#cortex <- readVcf(.testfile("cortex-1.0.5.14.vcf"), "")
#crest <- readVcf(.testfile("crest.vcf"), "")
delly <- readVcf(.testfile("delly-0.6.8.vcf"), "")
#gasv <- readVcf(.testfile("gasv-20140228.vcf"), "")
#gridss <- readVcf(.testfile("gridss-0.10.0.vcf"), "")
#lumpy <- readVcf(.testfile("lumpy-0.2.11.vcf"), "")
pindel <- readVcf(.testfile("pindel-0.2.5b6.vcf"), "")
#socrates <- readVcf(.testfile("socrates-1.13.vcf"), "")

# callers <- list()
# callers[["breakdancer"]] <- breakdancer
# callers[["cortex"]] <- cortex
# callers[["crest"]] <- crest
# callers[["delly"]] <- delly
# callers[["gasv"]] <- gasv
# callers[["gridss"]] <- gridss
# callers[["lumpy"]] <- lumpy
# callers[["pindel"]] <- pindel
# callers[["socrates"]] <- socrates


# unpaired breakend
test_that("mateBreakend fails if missing mate", {
    expect_error(mateBreakend(svToBreakendGRanges(simple)[1,]))
})

test_that("svToBreakendGRanges convert to breakend pairs", {
    gr <- svToBreakendGRanges(simple)
    pairId <- c("INS", "DEL", "SYMINS", "SYMDEL")
    expect_equal(mateBreakend(gr[paste0(pairId, "_bp1"),], gr[paste0(pairId, "_bp2"),]))
    #expect_equal(mateBreakend(gr[c("BNDFB", "BNDBF"),], gr[c("BNDBF", "BNDFB"),]))
})

test_that("svToBreakendGRanges non-symbolic alleles", {
    gr <- svToBreakendGRanges(simple[c("INS", "DEL"),])
    expect_equal(4, length(gr))
})