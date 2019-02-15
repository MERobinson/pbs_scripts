#!/user/bin/env bash
set -o pipefail
set -o nounset

# default args
workdir=$PWD
outdir=''
check='on'

# help message
help_message="
Collect insert size metrics with Picard CollectInsertSizeMetrics

usage:
    $(basename "$0") [-options] -b <BAM>
required arguments:
    -b|--bam : BAM file to analyse
optional arguments:
    -o|--outdir : output directory for metric files (default = pwd)
    -l|--logdir : output dir for log (default = --outdir)
    -n|--name : name prefix for output files (default = bam prefix)
    --check : whether to check input files [on|off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling

"

# parse arguments
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam=$2
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
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "ERROR: Undefined argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
	shift
done

# create outdir
if [[ -z ${outdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# check BAM arg and split
if [[ -z ${bam:-} ]]; then
	printf "\nERROR: --bam argument required\n"
	echo "$help_message"; exit 1
fi
if [[ ${check:-} = yes ]] && [[ ! -r ${workdir}/${bam} ]]; then
    printf "\nERROR: BAM file is not readable: %s/%s\n" $workdir $bam
    echo "$help_message"; exit 1
fi

# set name if not given
if [[ -z ${name:-} ]]; then
    name=$(basename ${bam})
    name=${name%%.*}
fi

# set commands
cmd=("picard CollectInsertSizeMetrics"
     "I=$(basename ${bam})"
     "O=${name}.isize_metrics.txt"
     "H=${name}.isize_histogram.pdf")

# run job
script=$(cat <<- EOS | qsub -N ${name}.isize -
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:ncpus=1:mem=8gb
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $logdir/${name}.$(basename $0 .sh).log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# load modules
		module load anaconda3/personal
		source activate picard
	
		# copy to scratch
		cp -L $workdir/${bam%%.*}* .

		# run cmd
		${cmd[@]}

		# copy metrics to outdir
		cp ${name}.isize_metrics.txt $workdir/$outdir/
		cp ${name}.isize_histogram.pdf $workdir/$outdir/
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		ls -lhAR
	EOS
)
echo "JOBID: $script"
exit 0
