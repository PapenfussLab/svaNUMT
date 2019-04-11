#' StructuralVariantAnnotation: a package for VCF
#'
#' plyranges is a dplyr like API to the Ranges/GenomicRanges infrastructure
#' in Bioconductor.
#'
#' StructuralVariantAnnotation contains useful helper functions for reading
#' and interpreting structural variants calls. The packages contains functions
#' for parsing VCFs from a number of popular caller as well as functions for
#' dealing with breakpoints involving two separate genomic loci. The package
#' supports encoding breakpoints as `GRanges` objects and BEDPE-formatted
#' data frame.
#'
#'    * Parse VCF objects with the `breakpointRanges()` function.
#'    * Find breakpoint overlaps with the `findBreakpointOverlaps()`
#'   and `countBreakpointOverlaps()` functions.
#'    * Generate BEDPE files for circos plot with `breakpointgr2bedpe()` function.
#'    * ...
#'
#' For more details on the features of StructuralVariantAnnotation, read the vignette:
#' `browseVignettes(package = "StructuralVariantAnnotation")`
#'
#' @docType package
#' @name StructuralVariantAnnotation
#' @import assertthat
#' @import BiocGenerics
#' @import VariantAnnotation
#' @import stringr
#' @import GenomicRanges
#' @import Biostrings
#' @import rtracklayer
#' @import dplyr
NULL
