
expect_equal(elementExtract(CharacterList(c("a"), c("b", "c"), c(), c("d", "e", "f"), "g")), c("a", "b", NA_character_, "d", "g"))
expect_equal(elementExtract(CharacterList(c("b", "c"))), c("b"))
expect_equal(elementExtract(CharacterList()), character(0))
expect_equal(elementExtract(CharacterList(c())), c(NA_character_))
expect_equal(elementExtract(c("aa", "bb", "aa")), c("aa", "bb", "aa"))
expect_equal(elementExtract(c("aa", "bb", "aa"), 3), c(NA_character_, NA_character_, NA_character_))

expect_equal(elementExtract(DNAStringSetList(c("A"), c("C", "G"), c())), c("A", "C", NA_character_))
expect_equal(elementExtract(DNAStringSet(c("A", "C", "N"))), c("A", "C", "N"))

expect_equal(elementExtract(IntegerList(c(), c(1, NA), c(2, 3), c(4, 5))), c(NA, 1, 2, 4))
expect_equal(elementExtract(IntegerList(c(), c(1, NA), c(2, 3), c(4, 5)), 2), c(NA, NA, 3, 5))
expect_equal(elementExtract(IntegerList()), integer(0))

expect_equal(elementExtract(NULL), NULL)

expect_equal(NULL %na% c(1,2,3), c(1,2,3))
expect_equal(c(1,NA,NA) %na% c(3,2,NA), c(1,2,NA))
expect_equal(c(1,2,3) %na% NULL, c(1,2,3))
expect_equal(c(1,NA) %na% c(1,2,3), c(1,2))
expect_equal(c(1,NA) %na% integer(0), c(1,NA))
# want to be able to perform operations on VCF columns that may not exist
expect_equal(numeric(0) %na% c(1,2), c(1,2))
expect_equal((NULL - c(2, NA)) %na% c(1,2), c(1,2))
