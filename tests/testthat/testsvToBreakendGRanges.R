
example <- readVcf(.testfile("vcf4.2.example.sv.vcf"), "")
simple <- readVcf(.testfile("simple.vcf"), "")
breakend <- readVcf(.testfile("breakend.vcf"), "")

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
test_that("partner fails if missing mate", {
    expect_error(partner(svToBreakendGRanges(breakend)[1,]))
})

test_that("svToBreakendGRanges convert to breakend pairs", {
    gr <- svToBreakendGRanges(simple)
    pairId <- c("INS", "DEL", "SYMINS", "SYMDEL")
    expect_true(all(paste0(pairId, "_bp1") %in% names(gr)))
    expect_true(all(paste0(pairId, "_bp2") %in% names(gr)))
    expect_equal(names(partner(gr))[names(gr) %in% paste0(pairId, "_bp1")], paste0(pairId, "_bp2"))
	expect_equal(names(partner(gr))[names(gr) %in% c("BNDFB", "BNDBF")], c("BNDBF", "BNDFB"))
})

test_that("svToBreakendGRanges non-symbolic alleles", {
    gr <- svToBreakendGRanges(simple[c("INS", "DEL"),])
    expect_equal(4, length(gr))
    
    gr <- svToBreakendGRanges(.testrecord("chr1	1	.	ATT	AGGA	.	.	"))
    expect_equal(start(gr), c(1, 4))
    expect_equal(gr$insSeq, c("GGA", "GGA"))
    expect_equal(gr$svLen, 1)
    
    gr <- svToBreakendGRanges(.testrecord("chr1	2	.	TTT	AGGA	.	.	"))
    expect_equal(start(gr), c(1, 4))
    expect_equal(gr$insSeq, c("AGGA", "AGGA"))
    expect_equal(gr$svLen, 1)
    
    gr <- svToBreakendGRanges(.testrecord("chr1	2	.	AGT	AGGA	.	.	"))
    expect_equal(start(gr), c(3, 5))
    expect_equal(gr$insSeq, c("GA", "GA"))
    expect_equal(gr$svLen, 1)
})

test_that("svToBreakendGRanges breakend", {
    gr <- svToBreakendGRanges(breakend)
    expect_equal("parid_b", gr["parid_a",]$partner)
    expect_equal("mateid_b", gr["mateid_a",]$partner)
    expect_equal(partner(gr[c("parid_a", "parid_b"),], gr[c("parid_b", "parid_a"),]))
    expect_warning(svToBreakendGRanges(breakend[c("mateid_a","mateid_b","multi_mateid", )]), "Ignoring additional mate breakends")
    expect_warning(svToBreakendGRanges(breakend[c("unpaired", )]), "Removing [0-9]+ unpaired breakend variants")
})