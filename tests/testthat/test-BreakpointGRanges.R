
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

"chr1	9595527	gridss43448o	A	A[chr1:9597585[ 1502.22	.	AS=1;ASCRP=3;ASCSR=25;ASQ=487.6149;ASRP=36;ASSR=25;BQ=38.23;BUM=2;BUMQ=38.227955;CIGAR=92M1X229N1X;CIPOS=-115,115;CIRPOS=-115,115;CQ=1477.48;EVENT=gridss43448;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448h;RAS=1;RASQ=643.5017;REF=58;REFPAIR=28;RP=15;RPQ=371.10327;SVTYPE=BND",
"chr1	9597585	gridss43448h	C	]chr1:9595527]C 1502.22	.	AS=1;ASCRP=3;ASCSR=25;ASQ=643.5017;ASRP=36;ASSR=25;BQ=18.71;BSC=1;BSCQ=18.707733;CIGAR=1X229N1X30M1D198M;CIPOS=-115,115;CIRPOS=-115,115;CQ=1477.48;EVENT=gridss43448;HOMLEN=230;HOMSEQ=TGGGAGGCTGAGGCAGGCAGATCACTTGAGGCCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCCTGTCTCTACTAAAAATACAGAAAAATTAGCCAGGCATGGTGGCACGTGCCTGTAATCCAGCTACTCGTGAGGCAGAGGCAGGAGAATTGCTTGAACCCAGGAGGTGGAGGTTGCAGTGAGCTGAGATCATGCCACTGCACTCCAGCCTGGGTGACAGAG;MATEID=gridss43448o;RAS=1;RASQ=487.6149;REF=22;REFPAIR=11;RP=15;RPQ=371.10327;SVTYPE=BND"


