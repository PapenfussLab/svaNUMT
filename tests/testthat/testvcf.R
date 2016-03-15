
simple <- readVcf(.testfile("simple.vcf"), "")
multipleAlleles <- readVcf(.testfile("multipleAltSVs.vcf"), "")

expect_false(.hasSingleAllelePerRecord(multipleAlleles))
expect_true(.hasSingleAllelePerRecord(expand(multipleAlleles)))

expect_equal(
    c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    isSymbolic(simple))

expect_equal(
    c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    isStructural(simple))

expect_equal(
    c(0, 0, 1, -1, 1, -2, NA, NA, NA, NA),
    svLen(simple))

expect_true(isStructural(.testrecord("chr1	1	.	ATT	AGGA	.	.	")))
expect_false(isStructural(.testrecord("chr1	1	.	ATT	NNN	.	.	")))
