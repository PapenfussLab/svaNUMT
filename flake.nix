{
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      with import nixpkgs { inherit system; };
      with rPackages; {
        defaultPackage = buildRPackage {
          name = "NUMTDetect";
          src = ./.;
          propagatedBuildInputs = [
            GenomicRanges
            rtracklayer
            VariantAnnotation
            StructuralVariantAnnotation
            BiocGenerics
            assertthat
            Biostrings
            stringr
            dplyr
            rlang
            GenomeInfoDb
            S4Vectors
            GenomicFeatures
          ];
        };

      });
}
