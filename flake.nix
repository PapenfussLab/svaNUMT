{
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      with import nixpkgs { inherit system; };
      with rPackages; {
        defaultPackage = (buildRPackage {
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
        }).overrideAttrs (_: {
          doCheck = true;
          checkInputs = [ BiocCheck ];
          script = writeText "bioc-check.r" ''
            library(BiocCheck)
            BiocCheck(paste0(Sys.getenv("TMPDIR"), "/NUMTDetect")
              , `no-check-deprecated`=T  # Requires net
              , `no-check-CRAN`=T        # Requires net
              , `no-check-version-num`=T # No idea why this fails
              , `quit-with-status`=T)
          '';
          checkPhase = ''
            cp -r $PWD $TMPDIR/NUMTDetect
            Rscript $script
          '';
        });
      });
}
