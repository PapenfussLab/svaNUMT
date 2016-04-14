expect_equal(.as.matrix.list(CharacterList(c("a"), c("b", "c"), c(), c("d", "e", "f"), "g")),
			 matrix(c("a", NA, NA,
			 		  "b", "c", NA,
			 		  NA, NA, NA,
			 		  "d", "e", "f",
			 		  "g", NA, NA), ncol=3, byrow=TRUE))
expect_equal(.as.matrix.list(IntegerList(c(1), c(2), c(3))),
	matrix(c(1,2,3), ncol=1, byrow=TRUE))
expect_equal(.as.matrix.list(IntegerList(c(1,4), c(2,5), c(3,6))),
	matrix(c(1,4,2,5,3,6), ncol=2, byrow=TRUE))
expect_equal(.as.matrix.list(IntegerList(c(), c(2,5), c())),
	matrix(c(NA, NA, 2, 5, NA, NA), ncol=2, byrow=TRUE))
expect_equal(.as.matrix.list(IntegerList(c(), c(), c())),
	matrix(integer(0), ncol=0, nrow=3))
expect_equal(.as.matrix.list(CharacterList(c(), c(), c())),
	matrix(character(0), ncol=0, nrow=3))


