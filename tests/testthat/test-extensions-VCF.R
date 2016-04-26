
example <- readVcf(.testfile("vcf4.2.example.sv.vcf"), "")
simple <- readVcf(.testfile("simple.vcf"), "")
breakend <- readVcf(.testfile("breakend.vcf"), "")
multipleAlleles <- readVcf(.testfile("multipleAltSVs.vcf"), "")


breakdancer <- readVcf(.testfile("breakdancer-1.4.5.vcf"), "")
#cortex <- readVcf(.testfile("cortex-1.0.5.14.vcf"), "")
#crest <- readVcf(.testfile("crest.vcf"), "")
delly <- readVcf(.testfile("delly-0.6.8.vcf"), "")
#gasv <- readVcf(.testfile("gasv-20140228.vcf"), "")
gridss <- readVcf(.testfile("gridss.vcf"), "")
#lumpy <- readVcf(.testfile("lumpy-0.2.11.vcf"), "")
pindel <- readVcf(.testfile("pindel-0.2.5b6.vcf"), "")
#socrates <- readVcf(.testfile("socrates-1.13.vcf"), "")

test_that("Delly TRA", {
	# https://groups.google.com/forum/#!msg/delly-users/6Mq2juBraRY/BjmMrBh3GAAJ
	# Sorry, I forgot to post this to the delly-users list:
	# For a translocation, you have 2 double strand breaks, one on chrA and one on chrB.
	# This creates 4 "dangling" ends, chrA_left, chrA_right, chrB_left, chrB_right.
	# For a translocation you can join chrA_left with chrB_left (3to3), chrA_left with chrB_right (3to5),
	# chrA_right with chrB_left (5to3) and chrA_right with chrB_right (5to5).
	# In fact for a typical reciprocal translocation in prostate cancer (where two chromosomes exchange their end)
	# Delly calls 2 translocations at the breakpoint, one 3to5 and one 5to3. But obviously not all translocations are reciprocal.
	# -Tobias
    gr <- breakpointRanges(.testrecord(c("chr10	2991435	TRA00000001	N	<TRA>	.	LowQual	CIEND=0,100;CIPOS=0,50;SVTYPE=TRA;CHR2=chr1;END=19357517;CT=3to5;INSLEN=0")))
    expect_equal(2, length(gr))
    expect_equal(c(2991435, 19357517), start(gr))
    expect_equal(c(2991485, 19357617), end(gr))
    expect_equal(c("+", "-"), as.character(strand(gr)))
    expect_equal(c("chr10", "chr1"), as.character(seqnames(gr)))
    gr <- breakpointRanges(delly)
})
test_that("pindel RPL", {
	gr <- breakpointRanges(pindel[5,])
	expect_equal(c(1029142, 1029183), start(gr))
	expect_equal(c("+", "-"), as.character(strand(gr)))
	expect_equal(c(51, 51), gr$insLen)
})
test_that("empty VCF", {
	expect_equal(0, length(breakpointRanges(.testrecord(c()))))
})

expect_false(.hasSingleAllelePerRecord(multipleAlleles))
expect_true(.hasSingleAllelePerRecord(expand(multipleAlleles)))

test_that("isSymbolic", {
	expect_equal(
	    c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	    isSymbolic(simple))
})
test_that("isStructural", {
	expect_equal(
	    c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
	    isStructural(simple))
	expect_true(isStructural(.testrecord("chr1	1	.	ATT	AGGA	.	.	")))
	expect_false(isStructural(.testrecord("chr1	1	.	ATT	NNN	.	.	")))
})
test_that(".svLen", {
	expect_equal(c(0, 0, 1, -1, 1, -2, NA, NA, NA, NA), .svLen(simple))
	expect_equal(.svLen(.testrecord("chr1	100	.	A	<DEL>	.	.	SVLEN=-1")), c(-1))
	expect_equal(.svLen(.testrecord("chr1	100	.	A	<DUP>	.	.	END=101")), c(1))
})



# unpaired breakend
test_that("partner fails if missing mate", {
	expect_error(partner(breakpointRanges(breakend)[1,]))
})

test_that("breakpointRanges convert to breakend pairs", {
	gr <- breakpointRanges(simple)
	pairId <- c("INS", "DEL", "SYMINS", "SYMDEL")
	expect_true(all(paste0(pairId, "_bp1") %in% names(gr)))
	expect_true(all(paste0(pairId, "_bp2") %in% names(gr)))
	expect_equal(names(partner(gr))[names(gr) %in% paste0(pairId, "_bp1")], paste0(pairId, "_bp2"))
	expect_equal(names(partner(gr))[names(gr) %in% c("BNDFB", "BNDBF")], c("BNDBF", "BNDFB"))
})

test_that("breakpointRanges creates placeholder names", {
	expect_named(breakpointRanges(.testrecord(c(
		"chr1	100	.	A	<DEL>	.	.	SVLEN=-1",
		"chr1	100	.	A	<DEL>	.	.	SVLEN=-1"))))
})

test_that("breakpointRanges non-symbolic alleles", {
	gr <- breakpointRanges(simple[c("INS", "DEL"),])
	expect_equal(4, length(gr))

	gr <- breakpointRanges(.testrecord("chr1	1	.	ATT	AGGA	.	.	"))
	expect_equal(start(gr), c(1, 4))
	expect_equal(gr$insSeq, c("GGA", "GGA"))
	expect_equal(gr$svLen, c(1, 1))

	gr <- breakpointRanges(.testrecord("chr1	2	.	TTT	AGGA	.	.	"))
	expect_equal(start(gr), c(1, 5))
	expect_equal(gr$insSeq, c("AGGA", "AGGA"))
	expect_equal(gr$svLen, c(1, 1))

	gr <- breakpointRanges(.testrecord("chr1	2	.	AGT	AGGA	.	.	"))
	expect_equal(start(gr), c(3, 5))
	expect_equal(gr$insSeq, c("GA", "GA"))
	expect_equal(gr$svLen, c(1, 1))
})

test_that("breakpointRanges intervals", {
	# Position assumed to the left aligned
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=0"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(100, 101))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=1"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(101, 102))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=2"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(102, 103))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMSEQ=AAAAAAAAAA"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(110, 111))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	CIPOS=-5,10"))
	expect_equal(start(gr), c(95, 96))
	expect_equal(end(gr), c(110, 111))

	# CIPOS over homology
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	CIPOS=-5,10;HOMLEN=50;HOMESEQ=A"))
	expect_equal(start(gr), c(95, 96))
	expect_equal(end(gr), c(110, 111))

	# HOMLEN over HOMSEQ
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	AT	.	.	HOMLEN=50;HOMESEQ=A"))
	expect_equal(start(gr), c(100, 101))
	expect_equal(end(gr), c(150, 151))
})

test_that("breakpointRanges simple indel", {
	gr <- breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	SVLEN=-1"))
	expect_equal(start(gr), c(100, 102))

	gr <- breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	END=101"))
	expect_equal(start(gr), c(100, 102))

	breakpointRanges(.testrecord("chr1	100	.	A	<DEL:WITH:SUBTYPE>	.	.	END=101"))
	breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	SVTYPE=DEL:WITH:SUBTYPE;END=101"))

	# warning about incompatable SVLEN and END fields
	#expect_warning(breakpointRanges(.testrecord("chr1	100	.	A	<DEL>	.	.	END=101;SVLEN=-10")), "SVLEN")
})

test_that("breakpointRanges breakend", {
	gr <- breakpointRanges(breakend)
	expect_equal("parid_b", gr["parid_a",]$partner)
	expect_equal("mateid_b", gr["mateid_a",]$partner)
	expect_equal(partner(gr[c("parid_a", "parid_b"),]), gr[c("parid_b", "parid_a"),])
	expect_warning(breakpointRanges(breakend[c("mateid_a", "mateid_b", "multi_mateid")]), "Ignoring additional mate breakends")
	expect_warning(breakpointRanges(breakend[c("unpaired")]), "Removing [0-9]+ unpaired breakend variants")
	expect_equal(breakpointRanges(.testrecord(c(
		"chr1	100	a	N	N[chr1:105[	.	.	SVTYPE=BND;CIPOS=0,1;PARID=b",
		"chr1	105	b	N	]chr1:100]N	.	.	SVTYPE=BND;CIPOS=0,1;PARID=a"
	)))$svLen, c(-4, -4))
	expect_equal(breakpointRanges(.testrecord(c(
		"chr1	100	a	N	NAAAA[chr1:101[	.	.	SVTYPE=BND;CIPOS=0,1;PARID=b",
		"chr1	101	b	N	]chr1:100]TTTTN	.	.	SVTYPE=BND;CIPOS=0,1;PARID=a"
	)))$svLen, c(4, 4))
})

test_that("breakpointRanges INV", {
	# VCF example
	gr <- breakpointRanges(.testrecord("chr1	321682	INV0	T	<INV>	6	PASS	SVTYPE=INV;END=421681"))
	expect_equal(4, length(gr))
	expect_equal(start(gr), c(321682, 421682, 321681, 421681))
	expect_equal(as.character(strand(gr)), c("-", "-", "+", "+"))
	expect_equal(names(gr), c("INV0_bp1", "INV0_bp2", "INV0_bp3", "INV0_bp4"))

	gr <- breakpointRanges(.testrecord("chr1	321682	INV0	T	<INV>	6	PASS	SVTYPE=INV;END=421681;CIPOS=-2,1;CIEND=-3,4"))
	expect_equal(4, length(gr))
	expect_equal(start(gr), c(321680, 421679, 321679, 421678))
	expect_equal(end(gr), c(321683, 421686, 321682, 421685))

	expect_error(breakpointRanges(.testrecord("chr1	321682	INV0	T	<INV>	6	PASS	SVTYPE=INV")))
})

test_that("breakpointRanges DUP", {
	# VCF example
	gr <- breakpointRanges(.testrecord(c(
		"chr1	12665100	.	A	<DUP>	14	PASS	SVTYPE=DUP;END=12686200;SVLEN=21100",
		"chr1	18665128	.	T	<DUP:TANDEM>	11	PASS	SVTYPE=DUP;END=18665204;SVLEN=76")))
	expect_equal(4, length(gr))
	expect_equal(start(gr), c(12665100, 18665128,
							  12686200, 18665204))
	expect_equal(as.character(strand(gr)), c("-", "-", "+", "+"))

	gr <- breakpointRanges(.testrecord("chr1	12665100	.	A	<DUP>	14	PASS	SVTYPE=DUP;END=12686200;SVLEN=21100;CIPOS=-2,1;CIEND=-3,4"))
	expect_equal(2, length(gr))
	expect_equal(start(gr), c(12665098, 12686197))
	expect_equal(end(gr), c(12665101, 12686204))

	expect_error(breakpointRanges(.testrecord("chr1	321682	.	T	<DUP>	.	.	SVTYPE=DUP")))
})
