
expect_equal(unpack(CharacterList(c("a"), c("b", "c"), c(), c("d", "e", "f"), "g")),
			 matrix(c("a", NA, NA,
			 		  "b", "c", NA,
			 		  NA, NA, NA,
			 		  "d", "e", "f",
			 		  "g", NA), ncol=3, byrow=TRUE))

