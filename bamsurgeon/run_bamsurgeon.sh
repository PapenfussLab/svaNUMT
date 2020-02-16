#!/bin/sh

# adds up to 100 SNPs to a ~770 kb region around the LARGE gene
# requires samtools/bcftools

REF=/home/dong.rn/home/StructuralVariantAnnotation/GRIPper/references-hs37d5-hs37d5.fa
BAM=~/home/Papenfuss_lab/projects/sv_benchmark_old/data.HG002/HG002_hs37d5_60x_1.sc.bam
VAR=/home/dong.rn/home/StructuralVariantAnnotation/bamsurgeon/simsv_region.csv
BAM_OUTPUT=/home/dong.rn/home/StructuralVariantAnnotation/bamsurgeon/HG002_sim_sv.bam
INS_FASTA=/home/dong.rn/home/StructuralVariantAnnotation/bamsurgeon/simGene_transcripts.fa

command -v addsv.py >/dev/null 2>&1 || { echo "addsv.py isn't installed" >&2; exit 65; }

if [ ! -e $REF ]
then
    echo "can't find reference .fasta: $REF, please supply a bwa-indexed .fasta"
    exit 65
fi

if [ ! -e $REF.bwt ]
then
    echo "can't find $REF.bwt: is $REF indexed with bwa?"
    exit 65
fi

addsv.py -p 6 -v $VAR -f $BAM -r $REF -o $BAM_OUTPUT --aligner mem --keepsecondary --seed 1234 --inslib $INS_FASTA

if [ $? -ne 0 ]
then
  echo "addsv.py failed."
  exit 65
# else
#   echo "sv added. need to sort bam"
#   echo "sorting output bam..."
#   samtools sort -T ../test_data/testregion_sv_mut.sorted.bam -o ../test_data/testregion_sv_mut.sorted.bam ../test_data/testregion_sv_mut.bam
#   mv ../test_data/testregion_sv_mut.sorted.bam ../test_data/testregion_sv_mut.bam
#
#   echo "indexing output bam..."
#   samtools index ../test_data/testregion_sv_mut.bam
# fi
