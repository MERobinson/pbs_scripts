#!/bin/bash

name=$1

jobid=$(cat <<- EOS | qsub -N $name.anno_vcf -
	#!/bin/bash
	#PBS -l walltime=24:00:00
	#PBS -l select=1:ncpus=1:mem=16gb
    #PBS -j oe
	#PBS -o $name.variant_anno.log
	#PBS -q med-bio

	module load java/jdk-8u66
	module load picard/2.6.0
	module load gatk/3.7
	module load samtools/1.2

	cp -rL /work/mrobinso/PM_cell_lines/variants/$name.hc.vcf* .
	cp -rL /work/mrobinso/PM_cell_lines/resources/dbsnp_150_common.170710.hg38.vcf.gz* .
	cp -rL /work/mrobinso/PM_cell_lines/resources/genome* .

	java -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
		-R genome.fa \
		-T VariantAnnotator \
		--dbsnp dbsnp_150_common.170710.hg38.vcf.gz \
		-o $name.hc.anno.vcf \
		-V $name.hc.vcf

	cp $name.hc.anno.vcf /work/mrobinso/PM_cell_lines/
	EOS
)
echo "JOBID: $jobid"
