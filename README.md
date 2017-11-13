# StructuralVariantAnnotation

StructuralVariantAnnotation contains useful helper
functions for dealing with structural variants in VCF format.
The packages contains functions for parsing VCFs from a number
of popular callers as well as functions for dealing with 
breakpoints involving two separate genomic loci encoded as
GRanges objects.

This package is under active development and should be considered as an unstable devel release.
No guarantees regarding API stability are made and function signatures could change between versions.
I am very open to implementing suggestions and enhancements requests for features that you would
find useful in your analysis so please raise them using the issues tab above.

## Installation

As the StructuralVariantAnnotation package is not yet released on CRAN or BioConductor, the simplest way to install/update is to use the devtools package:

```
#install.packages("devtools")
library(devtools)
install_github("PapenfussLab/StructuralVariantAnnotation")
```

## Documentation

Basic R documentation is included but, as this package is still a work in progress, many functions have only skeletal documentation.

The following are the core functions provided by this package:

- `breakpointRanges(vcf)` for converting the SVs in a VCF to a GRanges object containing one entry per breakend. Breakends on the "+" strand indicate a break immediately after the given position, and "-" indicates a break immediately before the given position. As an example, a "-" at the start of a chromosome, connected to a "+" at the end of the chromosome indicates that the chromosome is circular.
- `findBreakpointOverlaps(gr)` the breakpoint equivalent of findOverlaps()
- `countBreakpointOverlaps(gr)` the breakpoint equivalent of countOverlaps()
- `partner(gr)` for returning the breakpoint partners for each breakend
- `unpack(vcf)` for converting the VCF INFO columns to a data frame
- `bedpe2breakpointgr(file)` for processing BEDPE files

NB: partner() as it is current implement is quite brittle and is merely a GRanges column gr\$partner that requires the partner breakend to exist in the same GRanges object. Not only is directly user-visible, but any operation that subsets (eg extract only breakends at genes) without including or excluding both breakends will break the GRanges for breakpoint purposes (as will changing gr\$partner or names(gr)).



