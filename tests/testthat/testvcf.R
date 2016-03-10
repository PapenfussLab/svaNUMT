
simple <- expand(readVcf(.testfile("simple.vcf"), ""))

expect_equal(
    c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    isSymbolic(simple))

expect_equal(
    c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    isStructural(simple))

