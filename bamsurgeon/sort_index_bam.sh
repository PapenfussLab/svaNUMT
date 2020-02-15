#!/bin/sh
echo "sv added. need to sort bam"
echo "sorting output bam..."
samtools sort -T HG002_sim_tmp.sorted.bam -o HG002_sim_sv.sorted.bam -@ 20 HG002_sim_sv.bam
echo "indexing output bam..."
samtools index -@ 5 HG002_sim_sv.sorted.bam
