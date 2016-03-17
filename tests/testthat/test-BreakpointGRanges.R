
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
#pindel <- readVcf(.testfile("pindel-0.2.5b6.vcf"), "")
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
    expect_error(partner(breakpointRanges(breakend)[1,]))
})


expect_equal(IRanges::width(.constrict(breakpointRanges(breakend))), rep(1, length(breakpointRanges(breakend))))
expect_equal(start(.constrict(breakpointRanges(.testrecord(c(
	"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;CIPOS=0,1;PARID=b",
	"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;CIPOS=0,1;PARID=a",
	"chr1	100000	c	N	N]chr1:100100]	.	.	SVTYPE=BND;CIPOS=0,1;PARID=d",
	"chr1	100100	d	N	N]chr1:100000]	.	.	SVTYPE=BND;CIPOS=0,1;PARID=c"
	))))), c(100000, 100100, # not 100001, 100101
		100000, 100101))

hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

test_that("breakpointSequence", {
	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		# CTC>   <TGC
		"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a"
		))), hg19, anchoredBases=3), c("CTCTGC", "GCAGAG"))
	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		# CTC> AC  <TGC
		"chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a"
		))), hg19, anchoredBases=3), c("CTCACTGC", "GCAGTGAG"))
	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		"chr12	1000000	a	N	[chr12:2000000[N	.	.	SVTYPE=BND;PARID=b", #GGATA
		"chr12	2000000	b	N	[chr12:1000000[N	.	.	SVTYPE=BND;PARID=a"  #GAGAA
		))), hg19, anchoredBases=5), c("TATCCGAGAA", "TTCTCGGATA"))
	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		"chr12	1000000	a	N	[chr12:2000000[N	.	.	SVTYPE=BND;PARID=b;CIPOS=-115,115",
		"chr12	2000000	b	N	[chr12:1000000[N	.	.	SVTYPE=BND;PARID=a;CIPOS=-115,115"
		))), hg19, anchoredBases=5), c("TATCCGAGAA", "TTCTCGGATA"))

	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		"chr1	9595627	a	A	A[chr1:9597590[	.	.	MATEID=b;SVTYPE=BND",
		"chr1	9597590	b	C	]chr1:9595627]C	.	.	MATEID=a;SVTYPE=BND"
		))), hg19, anchoredBases=5)[1], "CTCCAAATCC")
})
test_that("referenceSequence", {
	expect_equal(referenceSequence(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a"
		))), hg19, anchoredBases=3), c("CTCACT", "GCATGG"))
	expect_equal(referenceSequence(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a"
		))), hg19, anchoredBases=3), c("CTCACT", "GCATGG"))
	expect_equal(referenceSequence(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a"
		))), hg19, anchoredBases=3, followingBases=5), c("CTCACTAA", "GCATGGCG"))
	expect_equal(referenceSequence(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b;CIPOS=-5,5",
		"chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a;CIPOS=-5,5"
		))), hg19, anchoredBases=3, followingBases=5), c("CTCACTAA", "GCATGGCG"))
})

test_that("referenceHomology", {
	# TODO: why doesn't this call position make sense? It should be to 9597590
	expect_gte(referenceHomology(breakpointRanges(.testrecord(c(
			"chr1	9595527	gridss43448o	A	A[chr1:9597585[	1502.22	.	CIPOS=-115,115;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448h;SVTYPE=BND",
			"chr1	9597585	gridss43448h	C	]chr1:9595527]C	1502.22	.	CIPOS=-115,115;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448o;SVTYPE=BND"
		))), hg19)$homlen[1], 230)

	expect_lt(referenceHomology(breakpointRanges(.testrecord(c(
		"chr12	1000000	a	A	A[chr12:2000000[	.	.	MATEID=b;SVTYPE=BND",
		"chr12	2000000	b	C	]chr12:1000000]C	.	.	MATEID=a;SVTYPE=BND"
	))), hg19)$homlen[1], 3)
})






