#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default args
workdir=$PWD
outdir=''
check='yes'

# help message
help_message="
Wrapper to run Picard's CollectAlignmentSummaryMetrics.

usage:
    $(basename "$0") [-options] -b <bam,files>
required arguments:
    -b|--bam : BAM file to analyse
    -F|--fasta : reference FASTA
optional arguments:
    -o|--outdir : output directory for metric files (default = pwd)
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output files (default = bam filename)
    --check : whether to check input files [yes,no] (default = yes)
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

# set logdir & create outdirs
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $workdir/$logdir
mkdir -p $workdir/$outdir

# check BAM arg and split
if [[ -z ${bam:-} ]]; then
	printf "\nERROR: --bam argument required\n"
	echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
fi
if [[ ${check:-} = yes ]]; then
    if [[ ! -r ${workdir}/${bam} ]]; then
        printf "\nERROR: BAM file is not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r ${workdir}/${fasta} ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
fi

# set basenames
bam_base=$(basename "$bam")
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%%.*}

# set name if not given
if [[ -z ${name:-} ]]; then
    name=${bam_base%%.*}
fi

# set commands
picard_command=("java -Xmx16G -jar /apps/picard/2.6.0/picard.jar"
                "CollectAlignmentSummaryMetrics R=$fasta_base"
                "I=$bam_base O=$name.alignment_summary_metrics")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# run job
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:ncpus=1:mem=16gb
		#PBS -j oe
		#PBS -q med-bio
		#PBS -N $name.AS_metrics
		#PBS -o $std_log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load java/jdk-8u66 &>> $out_log
		module load picard/2.6.0 &>> $out_log
		module load samtools/1.2 &>> $out_log
	
		# copy to scratch
		cp -rL $workdir/$fasta_prefix* . &>> $out_log
		cp $workdir/$bam* . &>> $out_log

		# run
		printf "\nCollecting alignment summary metrics:\n" >> $out_log
		${picard_command[@]} &>> $out_log

		# copy metrics to outdir
		cp $name.alignment_summary_metrics $workdir/$outdir/ &>> $out_log
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR
		cp $out_log $workdir/$logdir/
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
