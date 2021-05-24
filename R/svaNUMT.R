#' svaNUMT: a package for NUMT detection
#'
#' svaNUMT contains functions for detecting NUMT events from structural variant 
#' calls. svaNUMT contains functions for detecting NUMT events from structural 
#' variant calls. It takes structural variant calls in GRanges of breakend 
#' notation and identifies NUMTs by nuclear-mitochondrial breakend junctions. 
#' The main function reports candidate NUMTs if there is a pair of valid 
#' insertion sites found on the nuclear genome within a certain distance 
#' threshold. The candidate NUMTs are reported by events.
#'
#' For more details on the features of StructuralVariantAnnotation, read the vignette:
#' `browseVignettes(package = "svaNUMT")`
#'
#' @docType package
#' @name svaNUMT
#' @import BiocGenerics
#' @import VariantAnnotation
#' @import rtracklayer
#' @import Biostrings
#' @import GenomicRanges
#' @import StructuralVariantAnnotation
#' @import S4Vectors
#' @import GenomeInfoDb
#' @import GenomicFeatures
#' @importFrom dplyr %>%
#' @importFrom methods as is setMethod setGeneric
#' @importFrom rlang .data
NULL
