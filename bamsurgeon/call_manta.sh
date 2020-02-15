#!/bin/bash

BAM=~/home/Papenfuss_lab/projects/StructuralVariantAnnotation/bamsurgeon/HG002_sim_sv.sorted.bam
REF=~/home/Papenfuss_lab/projects/StructuralVariantAnnotation/GRIPper/references-hs37d5-hs37d5.fa
DIR=~/home/Papenfuss_lab/projects/StructuralVariantAnnotation/bamsurgeon

#ln -s $BAM input.bam
#ln -s $BAM.bai input.bam.bai
/nix/store/cdxyb4j6yvb7lqk7wipa19pj2rfxd6si-manta-1.6.0/bin/configManta.py \
	--bam $BAM \
	--referenceFasta $REF \
	--runDir $DIR && \
$DIR/runWorkflow.py -m local -j 6  --quiet #&& \
#$SCRIPT_DIR/manta_combine.sh $DIR/results > $SCRIPT_DIR/hgsvc/$BAM_NAME.vcf

#DATA_DIR=~/home/Papenfuss_lab/projects/hgsvc/high_cov_alignment
#DIR=~/home/Papenfuss_lab/projects/SVEnsemble/hgsvc
#REF=~/home/Papenfuss_lab/projects/hgsvc/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
#SCRIPT_DIR=~/home/Papenfuss_lab/projects/SVEnsemble/SVEnsemble
#ls -alh $DIR
#for BAM in $DATA_DIR/*.cram.bam ; do
#	echo $BAM
#	BAM_NAME="`echo $BAM | cut -c 68-130`"
#	echo $BAM_NAME
#	rm -rf $DIR; mkdir $DIR 2>/dev/null; cd $DIR
#	ln -s $BAM input.bam
#	ln -s $BAM.bai input.bam.bai
#	#ls -ahl .
#	/nix/store/cdxyb4j6yvb7lqk7wipa19pj2rfxd6si-manta-1.6.0/bin/configManta.py \
#		--bam input.bam \
#		--referenceFasta $REF \
#		--runDir $DIR && \
#	$DIR/runWorkflow.py -m local -j 6  --quiet && \
#	$SCRIPT_DIR/manta_combine.sh $DIR/results > $SCRIPT_DIR/hgsvc/$BAM_NAME.vcf
#done
