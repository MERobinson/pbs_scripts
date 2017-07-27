#!/bin/bash
module load bcftools
workdir=$PWD

for bcf in $workdir/variants/*.merged.delly.bcf; do

    # get file name and condition
    name=$(basename $bcf)
    name=${name%.*}
    cond=${name%%.*}

    # convert to vcf
    bcftools view $bcf > $name.vcf    

    # set jexl expression
    jexl="! vc.getGenotype(\"${cond}\").isHomRef() && vc.getGenotype(\"${cond}\").getGQ() > 20"

    # run script
    jobid=$(cat <<- EOS | qsub -N $name.sv_filter - 
	#!/bin/bash
    	#PBS -l walltime=05:00:00
	#PBS -l select=1:mem=30gb:ncpus=1
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/logs/$name.filter_sv.log

	module load java/jdk-8u66
	module load gatk/3.6
	module load bcftools

	cp $workdir/$name.vcf* .
	cp $workdir/resources/genome* .

	# filter to remove variants with hom ref and GQ < 20
	java -Xmx32G -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R genome.fa \
		--variant $name.vcf \
		-select "! vc.getGenotype(\"${cond}\").isHomRef() && vc.getGenotype(\"${cond}\").getGQ() > 20" \
		-o $name.filtered.vcf

	# copy back to workdir
	cp $name.filtered.vcf* $workdir/variants/
	ls -lhRA

	EOS
    )
done

