
example <- readVcf(.testfile("vcf4.2.example.sv.vcf"), "")
simple <- readVcf(.testfile("simple.vcf"), "")
breakend <- readVcf(.testfile("breakend.vcf"), "")

# unpaired breakend
test_that("partner fails if missing mate", {
    expect_error(partner(breakpointRanges(breakend)[1,]))
})
requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly=FALSE)
hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

test_that(".constrict", {
	expect_equal(IRanges::width(.constrict(breakpointRanges(breakend))), rep(1, length(breakpointRanges(breakend))))
	expect_equal(start(.constrict(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;CIPOS=0,1;PARID=b",
		"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;CIPOS=0,1;PARID=a",
		"chr1	100000	c	N	N]chr1:100100]	.	.	SVTYPE=BND;CIPOS=0,1;PARID=d",
		"chr1	100100	d	N	N]chr1:100000]	.	.	SVTYPE=BND;CIPOS=0,1;PARID=c"
		))))), c(100000, 100100, # not 100001, 100101
			100000, 100101))
	gr <- .constrict(breakpointRanges(.testrecord(c(
		"chrM	1	a	G	]chrM:16571]G	.	.	SVTYPE=BND;PARID=b;CIPOS=-10,-5",
		"chrM	16571	b	G	G[chrM:1[	.	.	SVTYPE=BND;PARID=a;CIPOS=5,10"
		))), hg19)
	expect_equal(start(gr), c(1, 16571))
})


test_that("findBreakpointOverlaps", {
	expect_equal(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
			"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a",
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
		breakpointRanges(.testrecord(c(
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c")))),
		data.frame(queryHits=c(3,4), subjectHits=c(1,2)))

	expect_equal(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
			"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a",
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
		breakpointRanges(.testrecord(c(
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
		maxgap=2000),
		data.frame(queryHits=c(1,2,3,4), subjectHits=c(1,2,1,2)))

	expect_equal(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	1	a	N	N[chr1:100100[	.	.	SVTYPE=BND;CIPOS=0,100000;PARID=b",
			"chr1	100100	b	N	]chr1:1]N	.	.	SVTYPE=BND;PARID=a"))),
		breakpointRanges(.testrecord(c(
			"chr1	100000	c	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100100	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c")))),
		data.frame(queryHits=c(1,2), subjectHits=c(1,2)))
})


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
	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		"chr1	1	a	N	N[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
		"chr1	1	b	N	]chr1:1]N	.	.	MATEID=a;SVTYPE=BND"
	))), hg19, anchoredBases=2), c("NNNN", "NNNN"))
	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		"chr1	1	a	N	N[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
		"chr1	1	b	N	]chr1:1]N	.	.	MATEID=a;SVTYPE=BND"
	))), hg19, 0, 0), c("", ""))
	expect_equal(breakpointSequence(breakpointRanges(.testrecord(c(
		"chr1	1	a	N	N[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
		"chr1	1	b	N	]chr1:1]N	.	.	MATEID=a;SVTYPE=BND"
	))), hg19, 0, 1), c("N", "N"))
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

	expect_equal(referenceSequence(breakpointRanges(.testrecord(c(
		"chrM	1	a	G	]chrM:16571]G	.	.	SVTYPE=BND;PARID=b",
		"chrM	16571	b	G	G[chrM:1[	.	.	SVTYPE=BND;PARID=a"
		))), hg19, anchoredBases=2, followingBases=3), c("TCNNN", "TGNNN"))
})

test_that("referenceHomology", {
	expect_gte(referenceHomology(breakpointRanges(.testrecord(c(
			"chr1	9595527	gridss43448o	A	A[chr1:9597585[	1502.22	.	CIPOS=-115,115;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448h;SVTYPE=BND",
			"chr1	9597585	gridss43448h	C	]chr1:9595527]C	1502.22	.	CIPOS=-115,115;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448o;SVTYPE=BND"
		))), hg19)$inexacthomlen[1], 230)

	expect_lt(referenceHomology(breakpointRanges(.testrecord(c(
		"chr12	1000000	a	A	A[chr12:2000000[	.	.	MATEID=b;SVTYPE=BND",
		"chr12	2000000	b	C	]chr12:1000000]C	.	.	MATEID=a;SVTYPE=BND"
	))), hg19)$inexacthomlen[1], 3)
	expect_equal(referenceHomology(breakpointRanges(.testrecord(c(
		"chr1	1	a	A	A[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
		"chr1	1	b	C	]chr1:1]C	.	.	MATEID=a;SVTYPE=BND"
	))), hg19, 5)$inexacthomlen, c(NA,NA))
	expect_true(is.na(referenceHomology(breakpointRanges(.testrecord(c(
		"chr1	1	a	A	A[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
		"chr1	1	b	C	]chr1:1]C	.	.	MATEID=a;SVTYPE=BND",
		"chr12	1000000	aa	A	A[chr12:2000000[	.	.	MATEID=bb;SVTYPE=BND",
		"chr12	2000000	bb	C	]chr12:1000000]C	.	.	MATEID=aa;SVTYPE=BND"
	))), hg19, 5)$inexacthomlen[1]))
})

test_that("blastHomology", {
	#Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/usr/local/bioinf/bin", sep=":"))
	#bh <- blastHomology(gr, hg19, "~/blastdb/16SMicrobial")

})



#bh <- blastHomology(gr[c("gridss14536o", "gridss14536h"),], hg19, "~/blastdb/16SMicrobial")








