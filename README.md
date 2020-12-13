# RTDetect
# StructuralVariantAnnotation

StructuralVariantAnnotation contains useful helper
functions for dealing with structural variants in VCF format.
The packages contains functions for parsing VCFs from a number
of popular callers as well as functions for dealing with 
breakpoints involving two separate genomic loci encoded as
GRanges objects.

## Installation


The StructuralVariantAnnotation package can be installed using BioConductor:

```
if (!requireNamespace("BiocManager", quietly=TRUE)) {
	install.packages("BiocManager")
}
BiocManager::install("StructuralVariantAnnotation")
```

