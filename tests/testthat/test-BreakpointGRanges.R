context('parsing functions')
example <- readVcf(.testfile("vcf4.2.example.sv.vcf"), "")
simple <- readVcf(.testfile("simple.vcf"), "")
breakend <- readVcf(.testfile("breakend.vcf"), "")

test_that("bedpe2breakpointgr creates unique ids", {
	gr <- pairs2breakpointgr(import(.testfile("unnamed.bedpe")))
	expect_equal(c("1", "3", "2", "4"), as.character(seqnames(gr)))
	expect_equal(c(19356, 1300148+1, 19427, 1302837+1), start(gr))
	expect_equal(c(19356, 1300151, 19427, 1302840), end(gr))
	expect_equal(c("+", "+", "-", "+"), as.character(strand(gr)))
	expect_equal(c(1, 41, 1, 41), gr$QUAL)
	expect_equal(c("bedpe1_1", "bedpe2_1", "bedpe1_2", "bedpe2_2"), names(gr))
})

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
	expect_equal(as.data.frame(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
			"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a",
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
		breakpointRanges(.testrecord(c(
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))))),
		data.frame(queryHits=c(3,4), subjectHits=c(1,2)))

	expect_equal(as.data.frame(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
			"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a",
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
		breakpointRanges(.testrecord(c(
			"chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
		maxgap=2000, sizemargin=NULL)),
		data.frame(queryHits=c(1,2,3,4), subjectHits=c(1,2,1,2)))

	expect_equal(as.data.frame(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	1	a	N	N[chr1:100100[	.	.	SVTYPE=BND;CIPOS=0,100000;PARID=b",
			"chr1	100100	b	N	]chr1:1]N	.	.	SVTYPE=BND;PARID=a"))),
		breakpointRanges(.testrecord(c(
			"chr1	100000	c	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100100	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))))),
		data.frame(queryHits=c(1,2), subjectHits=c(1,2)))
	expect_equal(as.data.frame(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	10000	a	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=-1000"))),
		breakpointRanges(.testrecord(c(
			"chr1	10000	c	N	N[chr1:11000[	.	.	SVTYPE=BND;PARID=d",
			"chr1	11000	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
		maxgap=1)),
		data.frame(queryHits=c(1,2), subjectHits=c(1,2)))
	expect_equal(as.data.frame(findBreakpointOverlaps(breakpointRanges(.testrecord(c(
			"chr1	10000	a	N	<DEL>	.	.	SVTYPE=DEL;SVLEN=-1000"))),
		breakpointRanges(.testrecord(c(
			"chr1	11000	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c",
			"chr1	10000	c	N	N[chr1:11000[	.	.	SVTYPE=BND;PARID=d"))),
		maxgap=1)),
		data.frame(queryHits=c(1,2), subjectHits=c(2,1)))
})
test_that("findBreakpointOverlaps_match_both_sides", {
	gr1 = GRanges(seqnames="1", ranges=IRanges(start=c(1, 100), width=1), strand="+", partner=c("2", "1"))
	names(gr1) = c("1", "2")
	gr2 = GRanges(seqnames="1", ranges=IRanges(start=c(100, 100), width=1),  strand="+", partner=c("2", "1"))
	names(gr2) = c("1", "2")
	expect_equal(nrow(as.data.frame(findBreakpointOverlaps(gr1, gr2, sizemargin=NULL))), 0)
})

test_that("findBreakpointOverlaps: delly vs truth", {
	grdelly <- breakpointRanges(.testrecord(c(
		"chr12	6905246	DEL00000292	N	<DEL>	.	PASS	PRECISE;CIEND=-52,52;CIPOS=-52,52;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.6.8;CHR2=chr12;END=100419669;CT=3to5;INSLEN=0;PE=28;MAPQ=34;SR=30;SRQ=1;CONSENSUS=GTGCATACATTTCAGTGACCCGTTTTAGAAACAGAATTAATATGGTGAATAGAGAAAGAAGAAATCAGTGACTTTGGCCAGGCACAGTAGCTCACATCTGTAATCCCAGCACTTTGGGAGGCTGAGACAGTTGGTTGCTTGAGCCCAGGAGT"
		)))
	grtruth <- breakpointRanges(.testrecord(c(
		"chr12	6905247	truth_8545_o	C	C[chr12:100419669[	.	.	EVENT=truth_8545_;MATEID=truth_8545_h;PARID=truth_8545_h;SVTYPE=BND",
		"chr12	100419669	truth_8545_h	G	]chr12:6905247]g	.	.	EVENT=truth_8545_;MATEID=truth_8545_o;PARID=truth_8545_o;SVTYPE=BND"
		)))
	hits <- findBreakpointOverlaps(grdelly, grtruth, maxgap=200, ignore.strand=TRUE)
	expect_equal(2, nrow(as.data.frame(hits)))
})

test_that("countBreakpointOverlaps", {
  expect_equal(countBreakpointOverlaps(breakpointRanges(.testrecord(c(
    "chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
    "chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a",
    "chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
    "chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
    breakpointRanges(.testrecord(c(
      "chr1	200000	a	N	N[chr1:200100[	.	.	SVTYPE=BND;PARID=b",
      "chr1	200100	b	N	]chr1:200000]N	.	.	SVTYPE=BND;PARID=a",
      "chr1	100000	c	N	N[chr1:100200[	.	.	SVTYPE=BND;PARID=d",
      "chr1	100200	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c")))),
    c(0,0,1,1))
})
test_that("countBreakpointOverlaps uniqueAllocation", {
  expect_equal(countBreakpointOverlaps(breakpointRanges(.testrecord(c(
    "chr1	100000	a	N	N[chr1:100100[	1	.	SVTYPE=BND;PARID=b",
    "chr1	100100	b	N	]chr1:100000]N	1	.	SVTYPE=BND;PARID=a",
    "chr1	100000	c	N	N[chr1:100100[	2	.	SVTYPE=BND;PARID=d",
    "chr1	100100	d	N	]chr1:100000]N	2	.	SVTYPE=BND;PARID=c"))),
    breakpointRanges(.testrecord(c(
      "chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
      "chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a",
      "chr1	100000	c	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=d",
      "chr1	100100	d	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=c"))),
    countOnlyBest=TRUE), c(0,0,2,2))
})

test_that("extractBreakpointSequence", {
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    # CTC>   <TGC
    "chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
    "chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a"
  ))), hg19, anchoredBases=3), c("CTCTGC", "GCAGAG"))
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    # CTC> AC  <TGC
    "chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
    "chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a"
  ))), hg19, anchoredBases=3), c("CTCACTGC", "GCAGTGAG"))
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    "chr12	1000000	a	N	[chr12:2000000[N	.	.	SVTYPE=BND;PARID=b", #GGATA
    "chr12	2000000	b	N	[chr12:1000000[N	.	.	SVTYPE=BND;PARID=a"  #GAGAA
  ))), hg19, anchoredBases=5), c("TATCCGAGAA", "TTCTCGGATA"))
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    "chr12	1000000	a	N	[chr12:2000000[N	.	.	SVTYPE=BND;PARID=b;CIPOS=-115,115",
    "chr12	2000000	b	N	[chr12:1000000[N	.	.	SVTYPE=BND;PARID=a;CIPOS=-115,115"
  ))), hg19, anchoredBases=5), c("TATCCGAGAA", "TTCTCGGATA"))
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    "chr1	9595627	a	A	A[chr1:9597590[	.	.	MATEID=b;SVTYPE=BND",
    "chr1	9597590	b	C	]chr1:9595627]C	.	.	MATEID=a;SVTYPE=BND"
  ))), hg19, anchoredBases=5)[1], "CTCCAAATCC")
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    "chr1	1	a	N	N[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
    "chr1	1	b	N	]chr1:1]N	.	.	MATEID=a;SVTYPE=BND"
  ))), hg19, anchoredBases=2), c("NNNN", "NNNN"))
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    "chr1	1	a	N	N[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
    "chr1	1	b	N	]chr1:1]N	.	.	MATEID=a;SVTYPE=BND"
  ))), hg19, 0, 0), c("", ""))
  expect_equal(extractBreakpointSequence(breakpointRanges(.testrecord(c(
    "chr1	1	a	N	N[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
    "chr1	1	b	N	]chr1:1]N	.	.	MATEID=a;SVTYPE=BND"
  ))), hg19, 0, 1), c("N", "N"))
})
test_that("extractReferenceSequence", {
  expect_equal(extractReferenceSequence(breakpointRanges(.testrecord(c(
    "chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
    "chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a"
  ))), hg19, anchoredBases=3), c("CTCACT", "GCATGG"))
  expect_equal(extractReferenceSequence(breakpointRanges(.testrecord(c(
    "chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
    "chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a"
  ))), hg19, anchoredBases=3), c("CTCACT", "GCATGG"))
  expect_equal(extractReferenceSequence(breakpointRanges(.testrecord(c(
    "chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
    "chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a"
  ))), hg19, anchoredBases=3, followingBases=5), c("CTCACTAA", "GCATGGCG"))
  expect_equal(extractReferenceSequence(breakpointRanges(.testrecord(c(
    "chr1	100000	a	N	NAC[chr1:100100[	.	.	SVTYPE=BND;PARID=b;CIPOS=-5,5",
    "chr1	100100	b	N	]chr1:100000]ACN	.	.	SVTYPE=BND;PARID=a;CIPOS=-5,5"
  ))), hg19, anchoredBases=3, followingBases=5), c("CTCACTAA", "GCATGGCG"))
  
  expect_equal(extractReferenceSequence(breakpointRanges(.testrecord(c(
    "chrM	1	a	G	]chrM:16571]G	.	.	SVTYPE=BND;PARID=b",
    "chrM	16571	b	G	G[chrM:1[	.	.	SVTYPE=BND;PARID=a"
  ))), hg19, anchoredBases=2, followingBases=3), c("TCNNN", "TGNNN"))
})

test_that("calculateReferenceHomology", {
  expect_gte(calculateReferenceHomology(breakpointRanges(.testrecord(c(
    "chr1	9595527	gridss43448o	A	A[chr1:9597585[	1502.22	.	CIPOS=-115,115;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448h;SVTYPE=BND",
    "chr1	9597585	gridss43448h	C	]chr1:9595527]C	1502.22	.	CIPOS=-115,115;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448o;SVTYPE=BND"
  ))), hg19)$inexacthomlen[1], 230)
  
  expect_lt(calculateReferenceHomology(breakpointRanges(.testrecord(c(
    "chr12	1000000	a	A	A[chr12:2000000[	.	.	MATEID=b;SVTYPE=BND",
    "chr12	2000000	b	C	]chr12:1000000]C	.	.	MATEID=a;SVTYPE=BND"
  ))), hg19)$inexacthomlen[1], 3)
  expect_equal(calculateReferenceHomology(breakpointRanges(.testrecord(c(
    "chr1	1	a	A	A[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
    "chr1	1	b	C	]chr1:1]C	.	.	MATEID=a;SVTYPE=BND"
  ))), hg19, 5)$inexacthomlen, c(NA,NA))
  expect_true(is.na(calculateReferenceHomology(breakpointRanges(.testrecord(c(
    "chr1	1	a	A	A[chr1:1[	.	.	MATEID=b;SVTYPE=BND",
    "chr1	1	b	C	]chr1:1]C	.	.	MATEID=a;SVTYPE=BND",
    "chr12	1000000	aa	A	A[chr12:2000000[	.	.	MATEID=bb;SVTYPE=BND",
    "chr12	2000000	bb	C	]chr12:1000000]C	.	.	MATEID=aa;SVTYPE=BND"
  ))), hg19, 5)$inexacthomlen[1]))
})
test_that("pairs_round_trip", {
  for (f in c("gridss.bedpe", "unnamed.bedpe")) {
    pairs = import(system.file("extdata", f, package = "StructuralVariantAnnotation"))
    gr = pairs2breakpointgr(pairs)
    pairs2 = breakpointgr2pairs(gr)
    expect_equal(seqnames(first(pairs)), seqnames(first(pairs2)))
    expect_equal(start(first(pairs)), start(first(pairs2)))
    expect_equal(strand(first(pairs)), strand(first(pairs2)))
    expect_equal(seqnames(second(pairs)), seqnames(second(pairs2)))
    expect_equal(start(second(pairs)), start(second(pairs2)))
    expect_equal(strand(second(pairs)), strand(second(pairs2)))
    expect_equal(mcols(pairs)$name, mcols(pairs2)$name)
    expect_equal(mcols(pairs)$score, mcols(pairs2)$score)
  }
})
test_that("single_breakends_valididity", {
  colo829_vcf = VariantAnnotation::readVcf(system.file("extdata", "COLO829T.purple.sv.ann.vcf.gz", package = "StructuralVariantAnnotation"))
  colo829_bpgr <- breakpointRanges(colo829_vcf)
  colo829_begr <- breakendRanges(colo829_vcf)
  colo829_gr <- sort(c(colo829_begr, colo829_bpgr))
  .assertValidBreakpointGRanges(colo829_gr, allowSingleBreakends=TRUE)
  expect_error(.assertValidBreakpointGRanges(colo829_gr, allowSingleBreakends=FALSE))
})

#test_that("calculateBlastHomology", {
#	Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/usr/local/bioinf/bin", sep=":"))
#	bh <- calculateBlastHomology(gr, hg19, "~/blastdb/16SMicrobial")
#
#})
# test_that("performance_test_partner", {
# 	n = 10000
# 	gr = GRanges(
# 		seqnames="1",
# 		ranges=IRanges(start=1:(2*n), width=1),
# 		partner=c(paste0(1:n, "o"), paste0(1:n, "h")))
# 	names(gr)=c(paste0(1:n, "h"), paste0(1:n, "o"))
# 	tictoc::tic(paste0("Start", n))
# 	pgr = partner(gr)
# 	tictoc::toc()
# })
test_that("simpleEventLength", {
	expect_equal(simpleEventLength(breakpointRanges(.testrecord(c(
			"chr1	100000	a	N	]chr1:100009]N	.	.	SVTYPE=BND;PARID=b",
			"chr1	100009	b	N	N[chr1:100000[	.	.	SVTYPE=BND;PARID=a",
			"chr1	100000	c	N	NNNNNNNNNNN[chr1:100001[	.	.	SVTYPE=BND;PARID=d",
			"chr1	100001	d	N	]chr1:100000]NNNNNNNNNNN	.	.	SVTYPE=BND;PARID=c")))),
		c(10, 10, 10, 10))
	expect_equal(simpleEventLength(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	N[chr1:100002[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100002	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a")))),
		c(-1, -1))
})
test_that("simpleEventType", {
	expect_equal(simpleEventType(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a"
	)))), c("DEL", "DEL"))
	expect_equal(simpleEventType(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	]chr1:100100]N	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	N[chr1:100000[	.	.	SVTYPE=BND;PARID=a"
	)))), c("DUP", "DUP"))
	expect_equal(simpleEventType(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	NNNNNNNNNNNNNNNNNNNNNNNN[chr1:100001[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100001	b	N	]chr1:100000]NNNNNNNNNNNNNNNNNNNNNNNN	.	.	SVTYPE=BND;PARID=a"
	)))), c("INS", "INS"))
	expect_equal(simpleEventType(breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	NNNNNNNNNNNNNNNNNNNNNNNN[chr2:100001[	.	.	SVTYPE=BND;PARID=b",
		"chr2	100001	b	N	]chr1:100000]NNNNNNNNNNNNNNNNNNNNNNNN	.	.	SVTYPE=BND;PARID=a"
	)))), c("CTX", "CTX"))
})
test_that("findInsDupOverlaps", {
	gr1 = breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	]chr1:100009]N	.	.	SVTYPE=BND;PARID=b",
		"chr1	100009	b	N	N[chr1:100000[	.	.	SVTYPE=BND;PARID=a",
		"chr1	100000	c	N	NNNNNNNNNNN[chr1:100001[	.	.	SVTYPE=BND;PARID=d",
		"chr1	100001	d	N	]chr1:100000]NNNNNNNNNNN	.	.	SVTYPE=BND;PARID=c")))
	gr2 = breakpointRanges(.testrecord(c(
		"chr1	100000	a	N	N[chr1:100100[	.	.	SVTYPE=BND;PARID=b",
		"chr1	100100	b	N	]chr1:100000]N	.	.	SVTYPE=BND;PARID=a",
		"chr1	100000	c	N	NNNNNNNNNNN[chr1:100001[	.	.	SVTYPE=BND;PARID=d",
		"chr1	100001	d	N	]chr1:100000]NNNNNNNNNNN	.	.	SVTYPE=BND;PARID=c")))
	hits12 = findInsDupOverlaps(gr1, gr2, maxgap=1)
	hits21 = findInsDupOverlaps(gr2, gr1, maxgap=1)
	expect_equal(queryHits(hits12), c(1,2))
	expect_equal(subjectHits(hits12), c(3,4))
	# symmetrical
	expect_equal(subjectHits(hits12), queryHits(hits21))
	expect_equal(subjectHits(hits21), queryHits(hits12))
})







