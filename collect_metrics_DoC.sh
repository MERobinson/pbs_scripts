#!/bin/bash

# default args
WORKDIR=$PWD
BAMDIR=bam
LOGDIR=logs
RESDIR=resources
OUTDIR=qc/metrics
FASTA=genome.fa

# help message
USAGE="$(basename "$0") [-irlfn] -b <bam,files>

purpose:
-- simple wrapper to run GATK's DepthOfCoverage
required arguments:
	-b | --bam : comma separated list of BAM files
optional arguments:
	-i | --indir  : input directory to find bam files (default = bam)
	-r | --resdir : resources directory containing reference FASTA (default = resources)
	-l | --logdir : output directory for log files (default = logs)
	-o | --outdir : output directory for metrics files (default = qc/metrics)
	-f | --fasta : name of link to whole genome FASTA (default = genome.fa)
	-n | --name : name prefix for output files (default = extracted from first bam file)
	--intervals : intervals file of regions to measure depth (default = none)
	-omitBaseOutput : argument passed to DepthOfCoverage
example:
	bash $(basename "$0") -b bam1.bam,bam2.bam -n test --intervals exons.intervals"

# parse arg
while [[ $# -gt 1 ]]; do
        key=$1
	case $key in
                -b|--bam)
		BAM_LIST=$2
		shift
		;;
                -i|--indir)
		BAMDIR=$2
		shift
		;;
		-r|--resdir)
		RESDIR=$2
		shift
		;;
		-l|--logdir)
		LOGDIR=$2
		shift
		;;
		-o|--outdir)
		OUTDIR=$2
		shift
		;;
		--intervals)
		INTERVALS=$2
		shift
		;;
                -n|--name)
                NAME=$2
                shift
                ;;
		-omitBaseOutput)
		OMIT_ARG="-omitBaseOutput"
		shift
		;;
                -h|--help)
                echo "$USAGE"
                exit 1
                ;;
                *)
                echo "Error: Illegal argument"
                echo "$USAGE"
                exit 1
                ;;
        esac
	shift
done

# check BAM arg and split
if [[ -z "$BAM_LIST" ]]; then
	echo "Error: bam list required"
	echo "$USAGE"; exit 1
else
	IFS=',' read -r -a BAM_ARRAY <<< "$BAM_LIST"
fi

# create outdir if needed
mkdir -p $WORKDIR/$LOGDIR
mkdir -p $WORKDIR/$OUTDIR

# if intervals file provided set arg string
if [[ ! -z $INTERVALS ]]; then
	L_STRING=$(basename "$INTERVALS")
	L_STRING="-L $L_STRING"
	L_CP_ARG="cp $WORKDIR/$RESDIR/$INTERVALS ."
fi

# set name if not given
if [[ -z "$NAME" ]]; then
	NAME=${BAM_ARRAY[0]%.*}
fi

# run job
JOBID=$(cat <<- EOS | qsub -N $NAME.DoC -
	#!/bin/bash
	#PBS -l walltime=40:00:00
	#PBS -l select=1:ncpus=1:mem=32gb
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $WORKDIR/$LOGDIR/$NAME.depth_of_coverage.log

	# load modules
	module load java/jdk-8u66
	module load picard/2.6.0
	module load gatk/3.6
	module load samtools/1.2
	
	# copy to scratch
	cp -rL $WORKDIR/$RESDIR/$FASTA_BASE* .
	for BAM in ${BAM_ARRAY[@]}; do
		cp $WORKDIR/$BAMDIR/'$BAM'* .
		echo '$BAM' > bam.list
	done
	$L_CP_ARG
	
	# run
	java -Xmx32G -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
		-T DepthOfCoverage \
		-R $FASTA \
		-o $NAME \
		-I bam.list \
		-ct 10 -ct 20 -ct 30 -ct 40 \
		$L_STRING $OMIT_ARG
	
	cp $NAME* $OUTDIR/

	ls -lhAR
	EOS
)
