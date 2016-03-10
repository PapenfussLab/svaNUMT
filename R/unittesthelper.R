#' testthat helper utility to locate files used
#' for package tests
.testfile <- function(filename, location="extdata") {
    file <- system.file(location, filename, package="StructuralVariantAnnotation")
    if (!file.exists(file)) {
        file <- file.path(getwd(), "inst", location)
    }
    assert_that(file.exists(file))
    return(file)
}