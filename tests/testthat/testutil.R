

expect_equal(.elementExtract(CharacterList(c("a"), c("b", "c"), c(), c("d", "e", "f"), "g")), c("a", "b", NA, "d", "g"))
expect_equal(.elementExtract(CharacterList(c("b", "c"))), c("b"))
expect_equal(.elementExtract(CharacterList()), c())
expect_equal(.elementExtract(CharacterList(c())), c(NA_character_))
expect_equal(.elementExtract(c("aa", "bb", "aa")), c("aa", "bb", "aa"))
