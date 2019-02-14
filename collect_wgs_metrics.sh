#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default args
workdir=$PWD
outdir=''
check='yes'

# help message
help_message="

Wrapper to call Picard's CollectWgsMetrics

usage:
    bash $(basename "$0") [-options] -b <BAM> -F <FASTA>
required arguments:
    -b|--bam : aligned BAM file
    -F|--fasta : whole genome FASTA file
optional arguments:
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output directory for log files [string] (default = --outdir)
    -n|--name : prefix for output files [string] (default = extracted from BAM)
    -I|--intervals : intervals file of regions to analyse [intervals list] (default = NULL)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list of PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments used for job scheduling/pipelines

"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -I|--intervals)
            intervals=$2
            shift
            ;;
		--check)
            check=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "ERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z "${bam:-}" ]]; then
    printf "\nERROR: No --bam argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: No --fasta argument provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$bam ]]; then
        printf "\nERROR: BAM file is not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    elif [[ ! -z ${intervals:-} ]] & [[ ! -r $workdir/${intervals:-} ]]; then
        printf "\nERROR: Intervals list is not readable: %s/%s\n" $workdir $intervals
        echo "$help_message"; exit 1
    fi
fi

# create outdir if needed
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# if intervals file provided set arg string
if [[ ! -z ${intervals:-} ]]; then
	intervals_arg=$(basename "$intervals")
	intervals_arg="INTERVALS=${intervals_arg}"
	intervals_cp="cp $workdir/$intervals ."
fi

# set name if not given
if [[ -z "${name:-}" ]]; then
	name=$(basename "$bam")
    name=${name%%.*}
fi

# get basenames/prefixes
fasta_base=$(basename $fasta)
fasta_prefix=${fasta%.*}
bam_base=$(basename $bam)

# run job
jobid=$(cat <<- EOS | qsub -N $name.collect_wgs_metrics -
		#!/bin/bash
		#PBS -l walltime=10:00:00
		#PBS -l select=1:mem=20gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.collect_wgs_metrics.log
		${depend:-}

		# load modules
		module load java/jdk-8u66
		module load picard/2.6.0
		module load samtools/1.2
	
		# copy to scratch
		cp -rL $workdir/$fasta_prefix* .
		cp $workdir/$bam* .
		${intervals_cp:-}
	
		# call picard
		java -Xmx18G -jar /apps/picard/2.6.0/picard.jar \
			CollectWgsMetrics \
			I=$bam_base \
			O=$name.collect_wgs_metrics.txt \
			R=$fasta_base \
			${intervals_arg:-}

		cp $name.collect_wgs_metrics.txt $workdir/$outdir/

		ls -lhAR
	EOS
)
echo "JOBID: $jobid"
