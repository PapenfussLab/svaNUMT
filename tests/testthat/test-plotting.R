vcf.vcf <- readVcf(.testfile("vcf4.2.example.sv.vcf"), "")
vcf.gr <- breakpointRanges(vcf.vcf)
gridss.vcf <- readVcf(.testfile("gridss.vcf"), "")
gridss.gr <- breakpointRanges(gridss.vcf)

test_that("converting non-BND VCF is successful",{
	vcf.bedpe <- breakpointgr2bedpe(vcf.gr)
	expect_equal(c("1","2827693","2827695","1","2827708","2827710"),
				 as.character(vcf.bedpe[1,seq(1,6)]))
	expect_equal(c('2','321625','321702','2','321877','321950'),
				 as.character(vcf.bedpe[2,seq(1,6)]))
	expect_equal(c('rs2376870_bp1'), as.character(vcf.bedpe$name[1]))
	expect_equal(c('2:321682_T/<DEL>_bp1'), as.character(vcf.bedpe$name[2]))
	expect_equal(c(NA,'6'), as.character(vcf.bedpe$score[seq(1,2)]))
	expect_equal(c('+','+'), as.character(vcf.bedpe$strand1[seq(1,2)]))
	expect_equal(c('-','-'), as.character(vcf.bedpe$strand2[seq(1,2)]))
})

test_that("converting BND VCF is successful",{
	gridss.bedpe <- breakpointgr2bedpe(gridss.gr)
	expect_equal(c("chr1","chr12"), as.character(gridss.bedpe$chrom1))
	expect_equal(c("18992157","84349"), as.character(gridss.bedpe$start1))
	expect_equal(c("18992158","84350"), as.character(gridss.bedpe$end1))
	expect_equal(c("chr12","chr12"), as.character(gridss.bedpe$chrom2))
	expect_equal(c("84963532","4886680"), as.character(gridss.bedpe$start2))
	expect_equal(c("84963533","4886681"), as.character(gridss.bedpe$end2))
	expect_equal(c("gridss2o","gridss39o"), as.character(gridss.bedpe$name))
	expect_equal(c("55","627.96"), as.character(gridss.bedpe$score))
	expect_equal(c("+","-"), as.character(gridss.bedpe$strand1))
	expect_equal(c("-","+"), as.character(gridss.bedpe$strand2))
})

test_that("converting fails if partner breakpoint not available",{
	expect_error(breakpointgr2bedpe(gridss.gr[1,]))
})
