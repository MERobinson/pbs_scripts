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
	-b | --bam : BAM file to analyse
optional arguments:
	-i | --indir  : input directory to find bam file (default = bam)
	-r | --resdir : resources directory containing reference FASTA (default = resources)
	-l | --logdir : output directory for log files (default = logs)
	-o | --outdir : output directory for metrics files (default = qc/metrics)
	-f | --fasta : name of link to whole genome FASTA (default = genome.fa)
	-n | --name : name prefix for output files (default = extracted from first bam file)
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
                -n|--name)
                NAME=$2
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
if [[ -z "$BAM" ]] | [[ ! -e "$WORKDIR/$BAMDIR/$BAM" ]]; then
	echo "ERROR: bam input required"
	echo "$USAGE"; exit 1
fi

# create outdir if needed
mkdir -p $WORKDIR/$LOGDIR
mkdir -p $WORKDIR/$OUTDIR

# set name if not given
if [[ -z "$NAME" ]]; then
	NAME=${BAM%%.*}
fi

# run job
JOBID=$(cat <<- EOS | qsub -N $NAME.DoC -
	#!/bin/bash
	#PBS -l walltime=20:00:00
	#PBS -l select=1:ncpus=1:mem=32gb
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $WORKDIR/$LOGDIR/$NAME.alignment_summary_metrics.log

	# load modules
	module load java/jdk-8u66
	module load picard/2.6.0
	module load samtools/1.2
	
	# copy to scratch
	cp -rL $WORKDIR/$RESDIR/$FASTA_BASE* .
	cp $WORKDIR/$BAMDIR/$BAM* .
	
	# run
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar \
		CollectAlignmentSummaryMetrics \
		R=$FASTA \
		I=$BAM \
		O=$NAME.alignment_summary_metrics
	
	cp $NAME.alignment_summary_metrics $OUTDIR/
	ls -lhAR
	EOS
)
